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
 * Enhance fault attributes and estimates fault strikes and dips, 
 * by optimal surface voting.
 *
 * @author Xinming Wu, Univerity of Texas at Austin.
 * @version 2017.07.16
 */
public class OptimalSurfaceVoter {

  /**
   * Constructs a parallel optimal surface voter.
   * @param ru windwow in fault normal direction.
   * @param rv windwow in fault strike direction.
   * @param rw windwow in fault dip direction.
   */
  public OptimalSurfaceVoter(int ru, int rv, int rw) {
    _ru = ru;
    _rv = rv;
    _rw = rw;
    _lmin = -ru;
    _lmax =  ru;
    _nl = 1+_lmax-_lmin;
    updateShiftRanges();
    _rgf = new RecursiveGaussianFilter(max(rv,rw));
  }

  /**
   * Sets bound on fault surface slopes in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
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
  public void setSurfaceSmoothing(double usmooth1, double usmooth2) {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    updateSmoothingFilters();
  }

  /**
   * Computes a voting score map from an input fault attribute image 
   * and approximately estimated fault strikes and dips.
   * An optimal surface voting method to enhance a fault attribute 
   * image so that the noisy features (unrelated to faults) are suppressed 
   * while the fault features are cleaner and more continuous. 
   * In this method, we first automatically pick seed points from 
   * the input attribute image and use these seeds as control points 
   * to compute optimal surface patches that pass through the control 
   * points and follow globally maximum fault attribute values. 
   * We then consider all the computed surfaces as voters and define 
   * voting scores for each voter by using fault attribute values 
   * smoothed along the surface voter. We further accumulate voting 
   * scores of all the voters to compute a voting score map as a 
   * new fault attribute image, where fault features (with high scores) 
   * are much cleaner, sharper, and more continuous than those in 
   * the input attribute image.
   * @param d minimum radius between the seed points
   * @param fm threshold of attribute value for selecting seed points
   * @param ft array of thinned fault attribute (et., planarity)
   * @param pt array of thinned fault strike
   * @param tt array or thinned fault dip
   * @return arrays of a voting score map and more accurately estimated 
   * fault strike and dips.
   */
  public float[][][][] applyVoting(int d, float fm,
    float[][][] ft, float[][][] pt, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final FaultCell[] seeds = pickSeeds(d,fm,ft,pt,tt);
    final int ns = seeds.length;
    final int nu = _nl;
    final int nv = _rv*2+1;
    final int nw = _rw*2+1;
    final int[] ct = new int[1];
    final float[][][] fs = smooth(ft);
    final float[][][] fe = new float[n3][n2][n1];
    final float[][][] vp = new float[n3][n2][n1];
    final float[][][] vt = new float[n3][n2][n1];
    final float[][][] vm = new float[n3][n2][n1];
    Stopwatch sw = new Stopwatch();
    sw.start();
    Parallel.loop(ns,new Parallel.LoopInt() {
    public void compute(int is) {
      ct[0] += 1;
      if(ct[0]%1000==0)
        System.out.println("done: "+ct[0]+"/"+ns);
      FaultCell cell = seeds[is];
      float[][] rws = new float[3][nw];
      float[][] rvs = new float[3][nv];
      float[][] rus = new float[3][nu];
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      float[] u = cell.getFaultNormal();
      float[] v = cell.getFaultDipVector();
      float[] w = cell.getFaultStrikeVector();
      updateVectorMap(_ru,u,rus[0],rus[1],rus[2]);
      updateVectorMap(_rv,v,rvs[0],rvs[1],rvs[2]);
      updateVectorMap(_rw,w,rws[0],rws[1],rws[2]);
      surfaceVoting(i1,i2,i3,u,v,w,rus,rvs,rws,fs,fe,vp,vt,vm);
    }});
    double timeUsed = sw.time();
    System.out.println("time used: "+timeUsed+" seconds");
    normalization(fe);
    return new float[][][][]{fe,vp,vt};
  }

  /**
   * Finds seed points from a thinned fault attribute image.
   * We select seed candidates from the thinned attribute image by 
   * collecting all the image samples with attribute values that are 
   * larger than some threshold fm. We finally check in the true seeds 
   * from all the candidate points in the order from the one with highest 
   * attribute value to the one with the lowest value. In checking in 
   * the seeds, we compute the distances between the current candidate 
   * to all the previously checked in seeds and check the current 
   * candidate as a true seed only if the minimum distance is larger 
   * than some predefined radius d.
   * @param d minimum radius between the seed points
   * @param fm threshold of attribute value for selecting seed points
   * @param ft array of thinned fault attribute (et., planarity)
   * @param pt array of thinned fault strike
   * @param tt array or thinned fault dip
   * @return array of seed points.
   */
  public FaultCell[] pickSeeds(
    int d, float fm, float[][][] ft, float[][][] pt, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final ArrayList<FaultCell> cs = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float fti = ft[i3][i2][i1];
      float pti = pt[i3][i2][i1];
      float tti = tt[i3][i2][i1];
      if(fti>fm) {
        FaultCell cell = new FaultCell(i1,i2,i3,fti,pti,tti);
        cs.add(cell);
      }
    }}}
    int np = cs.size();
    int[] is = new int[np];
    float[] fs = new float[np];
    for (int ip=0; ip<np; ++ip) {
      is[ip] = ip;
      fs[ip] = cs.get(ip).getFl();
    }
    quickIndexSort(fs,is);
    int[][][] mark = new int[n3][n2][n1];
    ArrayList<FaultCell> seeds = new ArrayList<FaultCell>();
    for (int ip=np-1; ip>=0; --ip) {
      FaultCell cell = cs.get(is[ip]);
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      int b1 = i1-d; b1=max(b1,0);
      int b2 = i2-d; b2=max(b2,0);
      int b3 = i3-d; b3=max(b3,0);
      int e1 = i1+d; e1=min(e1,n1-1);
      int e2 = i2+d; e2=min(e2,n2-1);
      int e3 = i3+d; e3=min(e3,n3-1);
      boolean ok = true;
      for (int k3=b3;k3<=e3;k3++) {
      for (int k2=b2;k2<=e2;k2++) {
      for (int k1=b1;k1<=e1;k1++) {
        if(mark[k3][k2][k1]==1) {
          ok=false;
          break;
        }
      }}}
      if(ok) {
        seeds.add(cell);
        mark[i3][i2][i1] = 1;
      }
    }
    return seeds.toArray(new FaultCell[0]);
  }

  /**
   * Returns fault strike and dip angles for specified fault normal vectors.
   * The components w2 and w3 must not both be zero; that is, the
   * fault plane cannot be horizontal.
   * @param u1 array 1st component of fault normal vector.
   * @param u2 array 2nd component of fault normal vector.
   * @param u3 array 3rd component of fault normal vector.
   * @return fault strike and dip angles, in degrees.
   */
  public float[][][][] strikeAndDipFromNormal(
    float[][][] u1, float[][][] u2, float[][][] u3) {
    int n3 = u1.length;
    int n2 = u1[0].length;
    int n1 = u1[0][0].length;
    float[][][] fp = new float[n3][n2][n1];
    float[][][] ft = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float u1i = u1[i3][i2][i1];
        float u2i = u2[i3][i2][i1];
        float u3i = u3[i3][i2][i1];
        if(u1i>0f) {
          u1i = -u1i;
          u2i = -u2i;
          u3i = -u3i;
        }
        ft[i3][i2][i1] = toDegrees(acos(-u1i));
        fp[i3][i2][i1] = range360(toDegrees(atan2(-u3i,u2i)));
      }}
    }});
    return new float[][][][]{fp,ft};
  }

  public void getFaultValues(FaultSkin[] skins, byte[][][] fs) {
    byte k = 0;
    for(FaultSkin skin:skins) {
      k++;
      for(FaultCell cell:skin) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        fs[i3][i2][i1] = k;
    }}
  }

  /**
   * Thins fault images to include only ridges in fault likelihoods.
   * After thinning, may be only one voxel wide. Thinned fault strikes and
   * dips are set to zero where thinned fault likelihoods are zero.
   * @param flpt array {fl,fp,ft} of fault likelihoods, strikes, and dips.
   * @return array {fl,fp,ft} of thinned fault likelihoods, strikes, and dips.
   */
  public static float[][][][] thin(float[][][][] flpt) {
    float[][][] f = flpt[0];
    float[][][] p = flpt[1];
    float[][][] t = flpt[2];
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    f = copy(f);
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1.0);
    rgf.applyX0X(f,f);
    rgf.applyXX0(f,f);
    float[][][] ff = new float[n3][n2][n1];
    float[][][] pp = new float[n3][n2][n1];
    float[][][] tt = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      int i3m = max(i3-1,0);
      int i3p = min(i3+1,n3-1);
      for (int i2=0; i2<n2; ++i2) {
        int i2m = max(i2-1,0);
        int i2p = min(i2+1,n2-1);
        float[] fmm = f[i3m][i2m];
        float[] fm0 = f[i3m][i2 ];
        float[] fmp = f[i3m][i2p];
        float[] f0m = f[i3 ][i2m];
        float[] f00 = f[i3 ][i2 ];
        float[] f0p = f[i3 ][i2p];
        float[] fpm = f[i3p][i2m];
        float[] fp0 = f[i3p][i2 ];
        float[] fpp = f[i3p][i2p];
        float[] p00 = p[i3 ][i2 ];
        float[] t00 = t[i3 ][i2 ];
        for (int i1=0; i1<n1; ++i1) {
          float f000 = f00[i1];
          float p000 = p00[i1];
          float t000 = t00[i1];
          if ((                p000<= 22.5f && f0m[i1]<f000 && f0p[i1]<f000) ||
              ( 22.5f<=p000 && p000<= 67.5f && fpm[i1]<f000 && fmp[i1]<f000) ||
              ( 67.5f<=p000 && p000<=112.5f && fp0[i1]<f000 && fm0[i1]<f000) ||
              (112.5f<=p000 && p000<=157.5f && fpp[i1]<f000 && fmm[i1]<f000) ||
              (157.5f<=p000 && p000<=202.5f && f0p[i1]<f000 && f0m[i1]<f000) ||
              (202.5f<=p000 && p000<=247.5f && fmp[i1]<f000 && fpm[i1]<f000) ||
              (247.5f<=p000 && p000<=292.5f && fm0[i1]<f000 && fp0[i1]<f000) ||
              (292.5f<=p000 && p000<=337.5f && fmm[i1]<f000 && fpp[i1]<f000) ||
              (337.5f<=p000                 && f0m[i1]<f000 && f0p[i1]<f000)) {
            ff[i3][i2][i1] = f000;
            pp[i3][i2][i1] = p000;
            tt[i3][i2][i1] = t000;
            if(p000>180) p000=360-p000;
            if(p000>60&&p000<120) {ff[i3m][i2 ][i1] = f000;}
            if(p000<30&&p000>150) {ff[i3 ][i2m][i1] = f000;}
          } else {
            pp[i3][i2][i1] = NO_STRIKE;
            tt[i3][i2][i1] = NO_DIP;
          }
        }
      }
    }
    float[][][][] flptn = new float[][][][]{ff,pp,tt};
    return flptn;
  }

  /**
   * Picks an optimal surface within a small cube centered at a seed point, 
   * computes voting scores on the surface, and collects the voting scores 
   * in a voting score map.
   */
  public void surfaceVoting(
    int c1, int c2, int c3, float[] u, float[] v, float[] w,
    float[][] dus, float[][] dvs, float[][] dws, float[][][] fx, 
    float[][][] fe, float[][][] vp, float[][][] vt, float[][][] vm) {
    int nu = dus[0].length;
    int nv = dvs[0].length;
    int nw = dws[0].length;
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    // get samples in uvw box
    float[][][] fs = fillfloat(1f,nu,nv,nw);
    samplesInUvwBox(c1,c2,c3,dus,dvs,dws,fs,fx);
    // find the optimal fault surface in the uvw box
    float[][] sfu = findSurface(fs);
    ArrayList<Integer> k1s = new ArrayList<Integer>();
    ArrayList<Integer> k2s = new ArrayList<Integer>();
    ArrayList<Integer> k3s = new ArrayList<Integer>();
    float fa = 0.0f;
    for (int kw=0; kw<nw; ++kw) {
      float dw1 = dws[0][kw]+c1;
      float dw2 = dws[1][kw]+c2;
      float dw3 = dws[2][kw]+c3;
      for (int kv=0; kv<nv; ++kv) {
        float iu = sfu[kw][kv];
        int i1 = round(iu*u[0]+dvs[0][kv]+dw1);
        int i2 = round(iu*u[1]+dvs[1][kv]+dw2);
        int i3 = round(iu*u[2]+dvs[2][kv]+dw3);
        boolean inbox = true;
        if(i1< 0||i1>=n1  ) inbox = false;
        if(i2<=0||i2>=n2-1) inbox = false;
        if(i3<=0||i3>=n3-1) inbox = false;
        if(inbox) {
          k1s.add(i1);
          k2s.add(i2);
          k3s.add(i3);
          fa += fx[i3][i2][i1];
        }
      }
    }
    int np = k1s.size();
    fa /= np;
    boolean alignX2 = false;
    float[] pt = surfaceStrikeAndDip(u,v,w,sfu);
    if(abs(u[2])>abs(u[1])) alignX2 = true;
    for (int ip=0; ip<np; ++ip) {
      int i1 = k1s.get(ip);
      int i2 = k2s.get(ip);
      int i3 = k3s.get(ip);
      fe[i3][i2][i1] += fa;
      if (alignX2) {
        fe[i3-1][i2][i1] += fa;
        fe[i3+1][i2][i1] += fa;
        if(fa>vm[i3-1][i2][i1]) {
          vm[i3-1][i2][i1] = fa;
          vp[i3-1][i2][i1] = pt[0];
          vt[i3-1][i2][i1] = pt[1];
        }
        if(fa>vm[i3+1][i2][i1]) {
          vm[i3+1][i2][i1] = fa;
          vp[i3+1][i2][i1] = pt[0];
          vt[i3+1][i2][i1] = pt[1];
        }
      } else {
        fe[i3][i2-1][i1] += fa;
        fe[i3][i2+1][i1] += fa;
        if(fa>vm[i3][i2-1][i1]) {
          vm[i3][i2-1][i1] = fa;
          vp[i3][i2-1][i1] = pt[0];
          vt[i3][i2-1][i1] = pt[1];
        }
        if(fa>vm[i3][i2+1][i1]) {
          vm[i3][i2+1][i1] = fa;
          vp[i3][i2+1][i1] = pt[0];
          vt[i3][i2+1][i1] = pt[1];
        }
      }
      if(fa>vm[i3][i2][i1]) {
        vm[i3][i2][i1] = fa;
        vp[i3][i2][i1] = pt[0];
        vt[i3][i2][i1] = pt[1];
      }

    }
  }

  public float[][][] byteToFloat(byte[][][] gx) {
    int n3 = gx.length;
    int n2 = gx[0].length;
    int n1 = gx[0][0].length;
    float[][][] gf = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      gf[i3][i2][i1] = (float)gx[i3][i2][i1];
    }}}
    return gf;
  }

  private float[] surfaceStrikeAndDip(
    float[] u, float[] v, float[] w, float[][] sfu) {
    int nw = sfu.length;
    int nv = sfu[0].length;
    float[][] sfs = new float[nw][nv];
    _rgf.apply00(sfu,sfs);
    float a1 = 1f;
    float a2 = -0.5f*(sfs[_rw][_rv+1]-sfs[_rw][_rv-1]);
    float a3 = -0.5f*(sfs[_rw+1][_rv]-sfs[_rw-1][_rv]);
    float as = 1f/sqrt(1+a2*a2+a3*a3);
    a1 *= as;
    a2 *= as;
    a3 *= as;
    float t1 = u[0]*a1+v[0]*a2+w[0]*a3;
    float t2 = u[1]*a1+v[1]*a2+w[1]*a3;
    float t3 = u[2]*a1+v[2]*a2+w[2]*a3;
    float ts = 1f/sqrt(t1*t1+t2*t2+t3*t3);
    t1 *= ts;
    t2 *= ts;
    t3 *= ts;
    if (t1>0) {
      t1 = -t1;
      t2 = -t2;
      t3 = -t3;
    }
    float ft = toDegrees(acos(-t1));
    float fp = range360(toDegrees(atan2(-t3,t2)));
    return new float[]{fp,ft};
  }

  //for display only
  public float[][][] surfacePicking(FaultCell cell, float[][][] fx) {
    final int nu = _nl;
    final int nv = _rv*2+1;
    final int nw = _rw*2+1;
    float[][] rws = new float[3][nw];
    float[][] rvs = new float[3][nv];
    float[][] rus = new float[3][nu];
    int i1 = cell.getI1();
    int i2 = cell.getI2();
    int i3 = cell.getI3();
    float[] u = cell.getFaultNormal();
    float[] v = cell.getFaultDipVector();
    float[] w = cell.getFaultStrikeVector();
    updateVectorMap(_ru,u,rus[0],rus[1],rus[2]);
    updateVectorMap(_rv,v,rvs[0],rvs[1],rvs[2]);
    updateVectorMap(_rw,w,rws[0],rws[1],rws[2]);
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    // get samples in uvw box
    float[][][] fs = fillfloat(1f,nu,nv,nw);
    samplesInUvwBox(i1,i2,i3,rus,rvs,rws,fs,fx);
    float[][] sf = findSurface(fs);
    float[][][] fb = fillfloat(1,n1,n2,n3);
    for (int kw=0; kw<nw; ++kw) {
      float dw1 = rws[0][kw]+i1;
      float dw2 = rws[1][kw]+i2;
      float dw3 = rws[2][kw]+i3;
      for (int kv=0; kv<nv; ++kv) {
        float iu = sf[kw][kv];
        int k1 = round(iu*u[0]+rvs[0][kv]+dw1);
        int k2 = round(iu*u[1]+rvs[1][kv]+dw2);
        int k3 = round(iu*u[2]+rvs[2][kv]+dw3);
        boolean inbox = true;
        if(k1< 0||k1>=n1  ) inbox = false;
        if(k2<=0||k2>=n2-1) inbox = false;
        if(k3<=0||k3>=n3-1) inbox = false;
        if (inbox) fb[k3][k2][k1] = -10f;
      }
    }
    return fb;
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
  public float[][] findSurface(float[][][] fx) {
    final int n2 = fx.length;
    final int n1 = fx[0].length;
    final int nl = fx[0][0].length;
    final float[][] u = new float[n2][n1];
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothFaultAttributes(fx,fx);
    float[][] d = new float[n1][nl];
    for (int i2=0; i2<n2; ++i2) {
      accumulateForward(fx[i2],d);
      backtrackReverse(d,fx[i2],uf[i2]);
    }
    smoothSurface(u,u);
    return u;
  }

 /**
   * Smooths the specified surface. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of surface to be smoothed.
   * @param us output array of smoothed surface.
   */
  public void smoothSurface(float[][] u, float[][] us) {
    if (_ref1!=null) {
      _ref1.apply1(u,us);
    } else {
      copy(u,us);
    }
    if (_ref2!=null)
      _ref2.apply2(us,us);
  }

    /**
   * Smooths (and normalizes) fault attributes.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of fault attributes.
   * @param es output array[n2][n1][nl] of smoothed fault attributes.
   */
  public void smoothFaultAttributes(float[][][] fx, float[][][] fs) {
    smoothFaultAttributes1(_bstrain1,fx,fs);
    smoothFaultAttributes2(_bstrain2,fs,fs);
  }

  /**
   * Accumulates fault attributes in forward direction.
   * @param e input array of fault attributes.
   * @param d output array of accumulated fault attributes.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }

  /**
   * Returns the optimal path found by backtracking in reverse.
   * @param d array of accumulated fault attributes.
   * @param e array of fault attributes.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }

    /**
   * Computes a path by backtracking in reverse direction.
   * @param d input array of accumulated fault attributes.
   * @param e input array of fault attributes.
   * @param u output array of the path.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_lmin,d,e,u);
  }

  private void normalization(final float[][][] fx) {
    final int n3 = fx.length;
    final int n2 = fx[0].length;
    final int n1 = fx[0][0].length;
    sub(fx,min(fx),fx);
    final float fmax = 1f/max(fx);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] fx3 = fx[i3];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float fxi = 1-fx3[i2][i1]*fmax;
        fxi *= fxi; //fxi^2
        fxi *= fxi; //fxi^4
        fxi *= fxi; //fxi^8
        fx3[i2][i1] = 1-fxi;
      }}
    }});
  }


  private void updateVectorMap(
    int r, float[] u, float[] ru1, float[] ru2, float[] ru3) {
    for (int i=1; i<=r; ++i) {
      int kp = r+i;
      int km = r-i;
      float iu1 = i*u[0];
      float iu2 = i*u[1];
      float iu3 = i*u[2];
      ru1[kp] =  iu1;
      ru2[kp] =  iu2;
      ru3[kp] =  iu3;
      ru1[km] = -iu1;
      ru2[km] = -iu2;
      ru3[km] = -iu3;
    }
  }

  private void samplesInUvwBox(
    int c1, int c2, int c3,
    float[][] dus, float[][] dvs, float[][] dws, 
    float[][][] fb, float[][][] fx) {
    int n3 = fx.length;
    int n2 = fx[0].length;
    int n1 = fx[0][0].length;
    int nw = fb.length;
    int nv = fb[0].length;
    for (int kw=0; kw<nw; kw++) {
      float dw1 = dws[0][kw]+c1;
      float dw2 = dws[1][kw]+c2;
      float dw3 = dws[2][kw]+c3;
      for (int kv=0; kv<nv; kv++) {
        float dv1 = dw1+dvs[0][kv];
        float dv2 = dw2+dvs[1][kv];
        float dv3 = dw3+dvs[2][kv];
        int um = _lmins[kw][kv]+_ru;
        int up = _lmaxs[kw][kv]+_ru;
        for (int ku=um; ku<=up; ku++) {
          int i1 = round(dv1+dus[0][ku]);
          int i2 = round(dv2+dus[1][ku]);
          int i3 = round(dv3+dus[2][ku]);
          i1 = min(max(i1,0),n1-1);
          i2 = min(max(i2,0),n2-1);
          i3 = min(max(i3,0),n3-1);
          fb[kw][kv][ku] = 1-fx[i3][i2][i1];
        }
      }
    }
  }


  private float[][][] smooth(float[][][] ft) {
    int n3 = ft.length;
    int n2 = ft[0].length;
    int n1 = ft[0][0].length;
    float[][][] fs = new float[n3][n2][n1];
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(1);
    rgf.apply000(ft,fs);
    fs = sub(fs,min(fs));
    return mul(fs,1f/max(fs));
  }
  /**
   * Finds the optimal path by backtracking in accumulated fault attributes.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param lmin minimum lag corresponding to lag index zero.
   * @param d input array[ni][nl] of accumulated fault attributes.
   * @param e input array[ni][nl] of fault attributes.
   * @param u output array[ni] of the picked path.
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
   * @param b strain parameter in 1st dimension.
   * @param e input array of fault attributes to be smooothed.
   * @param es output array of smoothed fault attributes.
   */
  private static void smoothFaultAttributes1(
    int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    for (int i2=0; i2<n2; ++i2)
      smoothFaultAttributes1(bf,ef[i2],esf[i2]);
  }

    /**
   * Smooths fault attributes in 1st dimension.
   * Does not normalize attributes after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of fault attributes to be smooothed.
   * @param es output array of smoothed fault attributes.
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
   * Smooths fault attributes in 2nd dimension.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of fault attributes to be smooothed.
   * @param es output array of smoothed fault attributes.
   */
  private static void smoothFaultAttributes2(
    int b, float[][][] e, float[][][] es) {
    int bf = b;
    int n2 = e.length;
    int n1 = e[0].length;
    int nl = e[0][0].length;
    float[][] e1  = new float[n2][nl];
    float[][] es1 = new float[n2][nl];
    float[][] ef1 = new float[n2][nl];
    float[][] er1 = new float[n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  e[i2][i1];
        es1[i2] = es[i2][i1];
      }
      accumulate( 1,bf,e1,ef1);
      accumulate(-1,bf,e1,er1);
      for (int i2=0; i2<n2; ++i2)
      for (int il=0; il<nl; ++il)
        es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
    }
  }

  /**
   * Non-linear accumulation of fault attributes.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of fault attributes.
   * @param d output array[ni][nl] of accumulated fault attributes.
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
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2*_bstrain2);
  }

  private void updateShiftRanges() {
    int nw = _rw*2+1;
    int nv = _rv*2+1;
    _lmins = new int[nw][nv];
    _lmaxs = new int[nw][nv];
    for (int iw=-_rw; iw<=_rw; ++iw) {
    for (int iv=-_rv; iv<=_rv; ++iv) {
      float wv = sqrt(iw*iw+iv*iv);
      if(wv>2) {
        _lmins[iw+_rw][iv+_rv] = max(-round(wv),_lmin);
        _lmaxs[iw+_rw][iv+_rv] = min( round(wv),_lmax);
      }
    }}
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _nl; // number of lags
  private int _lmin,_lmax; // min,max lags
  private int[][] _lmins,_lmaxs;
  private int _ru,_rv,_rw;
  private int _esmooth = 1; // number of nonlinear smoothings of attributes
  private int _bstrain1 = 4; // inverse of bound on slope in 1st dimension
  private int _bstrain2 = 4; // inverse of bound on slope in 2nd dimension
  private double _usmooth1 = 2.0; // extent of smoothing surface in 1st dim
  private double _usmooth2 = 2.0; // extent of smoothing surface in 2nd dim
  private RecursiveExponentialFilter _ref1; // for smoothing surface
  private RecursiveExponentialFilter _ref2; // for smoothing surface
  private RecursiveGaussianFilter _rgf;
  private static final float NO_STRIKE = -0.00001f;
  private static final float NO_DIP    = -0.00001f;

}
