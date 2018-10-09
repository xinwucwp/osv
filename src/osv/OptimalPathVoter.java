/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package osv;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Enhance fault attributes and estimates fault strikes, and dips, 
 * by optimal surface voting.
 *
 * @author Xinming Wu, Univerity of Texas at Austin.
 * @version 2017.07.18
 */
public class OptimalPathVoter {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public OptimalPathVoter(int ru, int rv) {
    _ru = ru;
    _rv = rv;
    _lmin = -ru;
    _lmax =  ru;
    _nl = 1+_lmax-_lmin;
    updateShiftRanges();
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  /**
   * Sets bound on fault surface slopes in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
  }

  /**
   * Sets the number of nonlinear smoothings of fault attributes.
   * The default number of smoothings is one.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setAttributeSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

    /**
   * Sets extents of smoothing filters used to smooth an extracted fault surface.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setPathSmoothing(double usmooth1) {
    _usmooth1 = usmooth1;
    updateSmoothingFilters();
  }

  public float[][][] applyVoting(int d, float fm,
    float[][] ft, float[][] pt) {
    final int n2 = ft.length;
    final int n1 = ft[0].length;
    final FaultCell2[] seeds = pickSeeds(d,fm,ft,pt);
    final int ns = seeds.length;
    final int nu = _nl;
    final int nv = _rv*2+1;
    final int[] ct = new int[1];
    final float[][] fs = copy(ft);
    final float[][] fe = new float[n2][n1];
    final float[][] fc = new float[n2][n1];
    final float[][] w1 = new float[n2][n1];
    final float[][] w2 = new float[n2][n1];
    Stopwatch sw = new Stopwatch();
    sw.start();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      ct[0] += 1;
      if(ct[0]%1000==0)
        System.out.println("done: "+ct[0]+"/"+ns);
      FaultCell2 cell = seeds[is];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      float[] u = cell.getFaultNormal();
      float[] v = cell.getFaultStrikeVector();
      float[][] rvs = new float[2][nv];
      float[][] rus = new float[2][nu];
      updateVectorMap(_ru,u,rus[0],rus[1]);
      updateVectorMap(_rv,v,rvs[0],rvs[1]);
      pathVoting(i1,i2,u,v,rus,rvs,fs,fc,fe,w1,w2);
    }});
    double timeUsed = sw.time();
    System.out.println("time used: "+timeUsed+" seconds");
    normalizationAndPower(1,fe);
    return new float[][][] {fe,w1,w2};
  }

  public float[][][][] applyVotingAndPathes(int d, float fm,
    float[][] ft, float[][] pt) {
    final int n2 = ft.length;
    final int n1 = ft[0].length;
    final FaultCell2[] seeds = pickSeeds(d,fm,ft,pt);
    final int ns = seeds.length;
    final int nu = _nl;
    final int nv = _rv*2+1;
    final int[] ct = new int[1];
    final float[][] fs = copy(ft);
    final float[][] fe = new float[n2][n1];
    final float[][] fc = new float[n2][n1];
    final float[][] w1 = new float[n2][n1];
    final float[][] w2 = new float[n2][n1];
    final float[][] x1 = new float[ns][n1];
    final float[][] x2 = new float[ns][n1];
    final float[][][] fes = new float[ns][n2][n1];
    Stopwatch sw = new Stopwatch();
    sw.start();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      ct[0] += 1;
      if(ct[0]%1000==0)
        System.out.println("done: "+ct[0]+"/"+ns);
      FaultCell2 cell = seeds[is];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      float[] u = cell.getFaultNormal();
      float[] v = cell.getFaultStrikeVector();
      float[][] rvs = new float[2][nv];
      float[][] rus = new float[2][nu];
      updateVectorMap(_ru,u,rus[0],rus[1]);
      updateVectorMap(_rv,v,rvs[0],rvs[1]);
      float[][] xs = pathVotingx(i1,i2,u,v,rus,rvs,fs,fc,fes[is],w1,w2);
      x1[is] = xs[0];
      x2[is] = xs[1];
    }});
    double timeUsed = sw.time();
    System.out.println("time used: "+timeUsed+" seconds");
    //normalizationAndPower(fe);
    float[][][] xs = new float[][][]{x1,x2};
    return new float[][][][] {fes,xs};
  }

  public float[][][] accumulate(float[][][] fes) {
    int ns = fes.length;
    int n2 = fes[0].length;
    int n1 = fes[0][0].length;
    float[][][] fas = new float[ns][n2][n1];
    fas[0] = fes[0];
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(1);
    for (int is=1; is<ns; ++is) {
      ref.apply(fes[is],fes[is]);
      add(fas[is-1],fes[is],fas[is]);
    }
    /*
    for (int is=0; is<ns; ++is) {
      normalizationAndPower(0,fas[is]);
      normalizationAndPower(0,fes[is]);
    }
    */
    return fas;
  }

  public float[][][] getSubPathes(
    int c1, int c2,
    FaultCell2[] seeds, float[][] x1, float[][] x2) {
    ArrayList<float[]> x1a = new ArrayList<float[]>();
    ArrayList<float[]> x2a = new ArrayList<float[]>();
    int ns = seeds.length;
    for (int is=0; is<ns; ++is) {
      FaultCell2 seed = seeds[is];
      int i1 = seed.getI1();
      int i2 = seed.getI2();
      if(i1==c1&&i2==c2) {
        x1a.add(x1[is]);
        x2a.add(x2[is]);
      }
    }
    float[][] x1s = x1a.toArray(new float[0][]);
    float[][] x2s = x2a.toArray(new float[0][]);
    return new float[][][]{x1s,x2s};
  }

  // for displaying only
  public float[][][] getPathes(
    float[][] ep, float[][] ft, FaultCell2[] seeds) {
    final int ns = seeds.length;
    final int nu = _nl;
    final int nv = _rv*2+1;
    final int[] ct = new int[1];
    final int n2 = ep.length;
    final int n1 = ep[0].length;
    final float[][] fs = copy(ft);//smooth(ft);
    final float[][] x1 = new float[ns][nv];
    final float[][] x2 = new float[ns][nv];
    Stopwatch sw = new Stopwatch();
    sw.start();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      ct[0] += 1;
      if(ct[0]%1000==0)
        System.out.println("done: "+ct[0]+"/"+ns);
      FaultCell2 cell = seeds[is];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      float[] u = cell.getFaultNormal();
      float[] v = cell.getFaultStrikeVector();
      float[][] rvs = new float[2][nv];
      float[][] rus = new float[2][nu];
      updateVectorMap(_ru,u,rus[0],rus[1]);
      updateVectorMap(_rv,v,rvs[0],rvs[1]);
      // get samples in uvw box
      float[][] ft = fillfloat(1f,nu,nv);
      samplesInUvBox(i1,i2,rus,rvs,fs,ft);
      // find the optimal fault surface in the uvw box
      float[] pu = findPath(ft);
      int jv = 0;
      for (int kv=0; kv<nv; ++kv) {
        float iu = pu[kv];
        float x1i = iu*u[0]+rvs[0][kv]+i1;
        float x2i = iu*u[1]+rvs[1][kv]+i2;
        if(x1i<2) continue;
        if(x2i<0) continue;
        if(x1i>n1-2) continue;
        if(x2i>n2-1) continue;
        x1[is][jv] = x1i;
        x2[is][jv] = x2i;
        jv++;
      }
      x1[is] = copy(jv,0,x1[is]);
      x2[is] = copy(jv,0,x2[is]);
    }});
    double timeUsed = sw.time();
    System.out.println("time used: "+timeUsed+" seconds");
    return new float[][][] {x1,x2};
  }

  public FaultCell2[] getSeeds(
    int c1, int c2, float[][] ft, float[][] pt) {
    FaultCell2[] seeds = new FaultCell2[1];
    seeds[0] = new FaultCell2(c1,c2,ft[c2][c1],pt[c2][c1]);
    return seeds;
  }


  public float[][] thin(float[][] f, float[][] w1, float[][] w2) {
    int n2 = f.length;
    final int n1 = f[0].length;
    final float[][] ff = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float d1 = w1[i2][i1];
        float d2 = w2[i2][i1];
        float x1p = i1+d1;
        float x2p = i2+d2;
        float x1m = i1-d1;
        float x2m = i2-d2;
        float fi = f[i2][i1];
        float fp = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1p,x2p);
        float fm = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1m,x2m);
        if(fp<fi&&fm<fi) {
          ff[i2][i1] = fi;
        }
      }
    }});
    return ff;
  }

  public float[][] seedToImage(
    float fmin,
    int n1, int n2, FaultCell2[] cells, float[][] ep) {
    float[][] fs = new float[n2][n1];
    for (FaultCell2 cell:cells) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      if(ep[i2][i1]>fmin)
        fs[i2][i1] = ep[i2][i1];
    }
    return fs;
  }

  public float[][] seedToPoint(FaultCell2[] cells) {
    int nc = cells.length;
    float[][] xs = new float[2][nc];
    for (int ic=0; ic<nc; ++ic) {
      xs[0][ic] = cells[ic].getI1();
      xs[1][ic] = cells[ic].getI2();
    }
    return xs;
  }


  public FaultCell2[] pickSeeds(
    int d, float fm, float[][] ft, float[][] pt) {
    final int n2 = ft.length;
    final int n1 = ft[0].length;
    final ArrayList<FaultCell2> cs = new ArrayList<FaultCell2>();
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float fti = ft[i2][i1];
      float pti = pt[i2][i1];
      if(fti>fm) {
        FaultCell2 cell = new FaultCell2(i1,i2,fti,pti);
        cs.add(cell);
      }
    }}
    int np = cs.size();
    int[] is = new int[np];
    float[] fs = new float[np];
    for (int ip=0; ip<np; ++ip) {
      is[ip] = ip;
      fs[ip] = cs.get(ip).getFl();
    }
    quickIndexSort(fs,is);
    int[][] mark = new int[n2][n1];
    ArrayList<FaultCell2> seeds = new ArrayList<FaultCell2>();
    for (int ip=np-1; ip>=0; --ip) {
      FaultCell2 cell = cs.get(is[ip]);
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      boolean ok = true;
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k2][k1]==1) {
          ok=false;
          break;
        }
      }}
      if(ok) {
        seeds.add(cell);
        mark[i2][i1] = 1;
      }
    }
    return seeds.toArray(new FaultCell2[0]);
  }

  public void pathVoting(
    int c1, int c2, float[] u, float[] v,
    float[][] dus, float[][] dvs,
    float[][] fx, float[][] fc, float[][] fe, float[][] w1, float[][] w2) {
    int nu = dus[0].length;
    int nv = dvs[0].length;
    int n2 = fx.length;
    int n1 = fx[0].length;
    // get samples in uvw box
    float[][] fs = fillfloat(1f,nu,nv);
    samplesInUvBox(c1,c2,dus,dvs,fx,fs);
    // find the optimal fault surface in the uvw box
    float[] pu = findPath(fs);
    float fa = 0.0f;
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    ArrayList<Float> p1s = new ArrayList<Float>();
    ArrayList<Float> p2s = new ArrayList<Float>();
    for (int kv=0; kv<nv; ++kv) {
      float p1i = -1f;
      float p2i =  0f;
      if(kv==0) {
        p2i= pu[kv+1]-pu[kv];
      } else {
        p2i = pu[kv]-pu[kv-1];
      }
      float iu = pu[kv];
      int i1 = round(iu*u[0]+dvs[0][kv]+c1);
      int i2 = round(iu*u[1]+dvs[1][kv]+c2);
      boolean inbox = true;
      if(i1<=0||i1>=n1-1) inbox = false;
      if(i2<=0||i2>=n2-1) inbox = false;
      if(inbox) {
        k1s.add(i1);
        k2s.add(i2);
        fa += fx[i2][i1];
        float psi = 1f/sqrt(p2i*p2i+1);
        p1s.add(p1i*psi);
        p2s.add(p2i*psi);
      }
    }
    int np = k1s.size();
    fa /= np;
    float u1 = u[0];
    float u2 = u[1];
    float v1 = v[0];
    float v2 = v[1];
    for (int ip=0; ip<np; ++ip) {
      int i1 = k1s.get(ip);
      int i2 = k2s.get(ip);
      float p1i = p1s.get(ip);
      float p2i = p2s.get(ip);
      fe[i2][i1] += fa;
      if(fa>fc[i2][i1]) {
        fc[i2][i1] = fa;
        w1[i2][i1] = u1*p1i+v1*p2i;
        w2[i2][i1] = u2*p1i+v2*p2i;
      }
    }
  }

  public float[][] pathVotingx(
    int c1, int c2, float[] u, float[] v,
    float[][] dus, float[][] dvs,
    float[][] fx, float[][] fc, float[][] fe, float[][] w1, float[][] w2) {
    int nu = dus[0].length;
    int nv = dvs[0].length;
    int n2 = fx.length;
    int n1 = fx[0].length;
    // get samples in uvw box
    float[][] fs = fillfloat(1f,nu,nv);
    samplesInUvBox(c1,c2,dus,dvs,fx,fs);
    // find the optimal fault surface in the uvw box
    float[] pu = findPath(fs);
    float fa = 0.0f;
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    ArrayList<Float> p1s = new ArrayList<Float>();
    ArrayList<Float> p2s = new ArrayList<Float>();
    int jv = 0;
    float[] x1 = new float[nv];
    float[] x2 = new float[nv];
    for (int kv=0; kv<nv; ++kv) {
      float p1i = -1f;
      float p2i =  0f;
      if(kv==0) {
        p2i= pu[kv+1]-pu[kv];
      } else {
        p2i = pu[kv]-pu[kv-1];
      }
      float iu = pu[kv];
      float x1i = iu*u[0]+dvs[0][kv]+c1;
      float x2i = iu*u[1]+dvs[1][kv]+c2;
      int i1 = round(x1i);
      int i2 = round(x2i);
      boolean inbox = true;
      if(i1<=0||i1>=n1-1) inbox = false;
      if(i2<=0||i2>=n2-1) inbox = false;
      if(inbox) {
        k1s.add(i1);
        k2s.add(i2);
        fa += fx[i2][i1];
        float psi = 1f/sqrt(p2i*p2i+1);
        p1s.add(p1i*psi);
        p2s.add(p2i*psi);
        x1[jv] = x1i;
        x2[jv] = x2i;
        jv++;
      }
    }
    x1 = copy(jv,0,x1);
    x2 = copy(jv,0,x2);
    int np = k1s.size();
    fa /= np;
    float u1 = u[0];
    float u2 = u[1];
    float v1 = v[0];
    float v2 = v[1];
    for (int ip=0; ip<np; ++ip) {
      int i1 = k1s.get(ip);
      int i2 = k2s.get(ip);
      float p1i = p1s.get(ip);
      float p2i = p2s.get(ip);
      fe[i2][i1] = fa;
      if(fa>fc[i2][i1]) {
        fc[i2][i1] = fa;
        w1[i2][i1] = u1*p1i+v1*p2i;
        w2[i2][i1] = u2*p1i+v2*p2i;
      }
    }
    return new float[][]{x1,x2};
  }



  /**
   * Returns angle in range [0,360] degrees.
   * @param phi angle, in degrees.
   * @return angle in range [0,360] degrees.
   */
  public static float range360(double phi) {
    while (phi<0.0)
      phi += 360.0;
    while (phi>=360.0)
      phi -= 360.0;
    return (float)phi;
  }

    /**
   * Returns angle in range [-180,180] degrees.
   * @param phi angle.
   * @return angle in range [-180,180] degrees.
   */
  public static float range180(double phi) {
    while (phi<-180.0)
      phi += 360.0;
    while (phi>180.0)
      phi -= 360.0;
    return (float)phi;
  }

  /**
   * Extract optimal fault surface from an input fault attribute image.
   * @param fx input array for the fault attribute image.
   */
  public float[] findPath(float[][] fx) {
    final int n1 = fx.length;
    final int nl = fx[0].length;
    final float[] u = new float[n1];
    for (int is=0; is<_esmooth; ++is)
      smoothFaultAttributes(fx,fx);
    float[][] d = new float[n1][nl];
    accumulateForward(fx,d);
    backtrackReverse(d,fx,u);
    smoothPath(u,u);
    return u;
  }


 /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothPath(float[] u, float[] us) {
    if (_ref1!=null) {
      _ref1.apply(u,us);
    } 
  }

    /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothFaultAttributes(float[][] fx, float[][] fs) {
    smoothFaultAttributes1(_bstrain1,fx,fs);
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  /**
   * Returns shifts found by backtracking in reverse.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }

    /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_lmin,d,e,u);
  }

  private void normalizationAndPower(float sigma, final float[][] fx) {
    if(sigma>0) {
      RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(sigma);
      RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
      ref.apply(fx,fx);
      //rgf.apply00(fx,fx);
    }
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    sub(fx,min(fx),fx);
    final float fmax = 1f/max(fx);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[] fx2 = fx[i2];
      for (int i1=0; i1<n1; ++i1) {
        float fxi = 1-fx2[i1]*fmax;
        fxi *= fxi; //fxi^2
        fxi *= fxi; //fxi^4
        //fxi *= fxi; //fxi^8
        fx2[i1] = 1-fxi;
      }
    }});
  }

  private void normalization(final float[][] fx) {
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply00(fx,fx);
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    sub(fx,min(fx),fx);
    final float fmax = 1f/max(fx);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[] fx2 = fx[i2];
      for (int i1=0; i1<n1; ++i1)
        fx2[i1] *= fmax;
    }});
  }


  private void updateVectorMap(
    int r, float[] u, float[] ru1, float[] ru2) {
    for (int i=1; i<=r; ++i) {
      int kp = r+i;
      int km = r-i;
      float iu1 = i*u[0];
      float iu2 = i*u[1];
      ru1[kp] =  iu1;
      ru2[kp] =  iu2;
      ru1[km] = -iu1;
      ru2[km] = -iu2;
    }
  }

  private void samplesInUvBox(
    int c1, int c2,
    float[][] dus, float[][] dvs,
    float[][] fx, float[][] fb) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int nv = fb.length;
    for (int kv=0; kv<nv; kv++) {
      float dv1 = dvs[0][kv]+c1;
      float dv2 = dvs[1][kv]+c2;
      int um = _lmins[kv]+_ru;
      int up = _lmaxs[kv]+_ru;
      for (int ku=um; ku<=up; ku++) {
        int i1 = round(dv1+dus[0][ku]);
        int i2 = round(dv2+dus[1][ku]);
        i1 = min(max(i1,0),n1-1);
        i2 = min(max(i2,0),n2-1);
        fb[kv][ku] = 1-fx[i2][i1];
      }
    }
  }

  private void samplesInUvBox(
    int c1, int c2,
    float[][] dus, float[][] dvs,
    float[][] fx, float[][] fb1, float[][] fb2) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    int rv = fb1.length;
    for (int kv=0; kv<rv; kv++) {
      float dv1 = dvs[0][kv+rv]+c1;
      float dv2 = dvs[1][kv+rv]+c2;
      int um = _lmins[kv+rv]+_ru;
      int up = _lmaxs[kv+rv]+_ru;
      for (int ku=um; ku<=up; ku++) {
        int i1 = round(dv1+dus[0][ku]);
        int i2 = round(dv2+dus[1][ku]);
        i1 = min(max(i1,0),n1-1);
        i2 = min(max(i2,0),n2-1);
        fb1[kv][ku] = 1-fx[i2][i1];
      }
    }

    for (int kv=rv-1; kv>=0; kv--) {
      float dv1 = dvs[0][kv]+c1;
      float dv2 = dvs[1][kv]+c2;
      int um = _lmins[kv]+_ru;
      int up = _lmaxs[kv]+_ru;
      for (int ku=um; ku<=up; ku++) {
        int i1 = round(dv1+dus[0][ku]);
        int i2 = round(dv2+dus[1][ku]);
        i1 = min(max(i1,0),n1-1);
        i2 = min(max(i2,0),n2-1);
        fb2[rv-1-kv][ku] = 1-fx[i2][i1];
      }
    }

  }


  private float[][] smooth(float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    float[][] fs = new float[n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(2);
    rgf.apply00(ft,fs);
    fs = sub(fs,min(fs));
    return mul(fs,1f/max(fs));
  }

  private float[][] smooth(float fm, float[][] ft) {
    int n2 = ft.length;
    int n1 = ft[0].length;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      if(ft[i2][i1]<fm) ft[i2][i1] = 0f;
    }}
    return smooth(ft);
  }

  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void backtrack(
    int dir, int b, int lmin, float[][] d, float[][] e, float[] u) 
  {
    float ob = 1.0f/b;
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?nim1:0;
    int is = (dir>0)?1:-1;
    int ii = ib;
    int il = max(0,min(nlm1,-lmin));
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]<dl) {
        il = jl;
        dl = d[ii][jl];
      }
    }
    u[ii] = il+lmin;
    while (ii!=ie) {
      int ji = max(0,min(nim1,ii+is));
      int jb = max(0,min(nim1,ii+is*b));
      int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
      int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
      float dm = d[jb][ilm1];
      float di = d[ji][il  ];
      float dp = d[jb][ilp1];
      for (int kb=ji; kb!=jb; kb+=is) {
        dm += e[kb][ilm1];
        dp += e[kb][ilp1];
      }
      dl = min3(dm,di,dp);
      if (dl!=di) {
        if (dl==dm) {
          il = ilm1;
        } else {
          il = ilp1;
        }
      }
      ii += is;
      u[ii] = il+lmin;
      if (il==ilm1 || il==ilp1) {
        float du = (u[ii]-u[ii-is])*ob;
        u[ii] = u[ii-is]+du;
        for (int kb=ji; kb!=jb; kb+=is) {
          ii += is;
          u[ii] = u[ii-is]+du;
        }
      }
    }
  }

  /**
   * Smooths fault attributes in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothFaultAttributes1(
    int b, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,e,ef);
    accumulate(-1,b,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1;
    int ie = (dir>0)?ni:-1;
    int is = (dir>0)?1:-1;
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb = max(0,min(nim1,ii-is*b));
      for (int il=0; il<nl; ++il) {
        int ilm1 = il-1; if (ilm1==-1) ilm1 = 0;
        int ilp1 = il+1; if (ilp1==nl) ilp1 = nlm1;
        float dm = d[jb][ilm1];
        float di = d[ji][il  ];
        float dp = d[jb][ilp1];
        for (int kb=ji; kb!=jb; kb-=is) {
          dm += e[kb][ilm1];
          dp += e[kb][ilp1];
        }
        d[ii][il] = min3(dm,di,dp)+e[ii][il];
      }
    }
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
  }

  private void updateShiftRanges() {
    int nv = _rv*2+1;
    _lmins = new int[nv];
    _lmaxs = new int[nv];
    for (int iv=-_rv; iv<=_rv; ++iv) {
      if(abs(iv)>2) {
        _lmins[iv+_rv] = max(-abs(iv),_lmin);
        _lmaxs[iv+_rv] = min( abs(iv),_lmax);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private int[] _lmins,_lmaxs;
  private int _ru,_rv;
  private int _esmooth = 1; // number of nonlinear smoothings of attributes
  private int _bstrain1 = 4; // inverse of bound on slope in 1st dimension
  private double _usmooth1 = 2.0; // extent of smoothing shifts in 1st dim
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private static final float NO_STRIKE = -0.00001f;
  private static final float NO_DIP    = -0.00001f;
  private SincInterpolator _si;

}
