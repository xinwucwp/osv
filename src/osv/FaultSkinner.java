/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package osv;

import java.util.*;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Constructs fault surfaces, which are represented as fault skins.
 * @author Xinming Wu, Univerity of Texas at Austin.
 * @version 2017.07.16
 */
public class FaultSkinner {

  public void setMaxDeltaStrike(float dp) {
    _dpmax = dp;
  }

  public void setGrowing(int dw, float an) {
    _dw = dw;
    _an = an;
  }

  public FaultCell[] findSeeds(
    int d, float fm, float[][][] ep, float[][][] ft, 
    float[][][] pt, float[][][] tt) {
    final int n3 = ft.length;
    final int n2 = ft[0].length;
    final int n1 = ft[0][0].length;
    final ArrayList<FaultCell> cs = new ArrayList<FaultCell>();
    for (int i3=0; i3<n3; i3++) {
    for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      float epi = ep[i3][i2][i1];
      float fti = ft[i3][i2][i1];
      float pti = pt[i3][i2][i1];
      float tti = tt[i3][i2][i1];
      if(fti>fm&&epi>0.8f) {
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

  public FaultSkin[] findSkins(
    float fmin, int kmin, FaultCell[] seeds, 
    float[][][] fv, float[][][] fp, float[][][] ft) {
    _fp = fp;
    _ft = ft;
    int nseed = seeds.length;
    int n3 = fv.length;
    int n2 = fv[0].length;
    int n1 = fv[0][0].length;
    int rw = max(n2,n3);
    int rv = rw;
    int ru = 150;
    // Empty list of skins.
    ArrayList<FaultSkin> skinList = new ArrayList<FaultSkin>();
    FaultCellGrid fcgs = new FaultCellGrid(n1,n2,n3);
    _fcgx = new FaultCellGrid(n1,n2,n3);
    fcgs.set(seeds);
    for (int kseed=0; kseed<nseed; kseed++) {
      while (kseed<nseed && seeds[kseed].used)
        ++kseed;
      if(kseed<nseed) {
        System.out.println("kseed="+kseed);
        FaultCell seed = seeds[kseed];
        FaultSkin skin = findSkin(fmin,ru,rv,rw,seed,fv);
        if(skin.size()>kmin) {
          skinList.add(skin);
          _fcgx.set(skin);
          for (FaultCell cell:skin)
            fcgs.setCellsInBox(cell,5,5,5);
        }
      }
    }
    return skinList.toArray(new FaultSkin[0]);
  }

  public FaultSkin findSkin(
    float fmin, int ru, int rv, int rw, FaultCell seed, 
    float[][][] fv) {
    int n3 = fv.length;
    int n2 = fv[0].length;
    int n1 = fv[0][0].length;
    int nu = ru*2+1;
    int nv = rv*2+1;
    int nw = rw*2+1;
    float o1 = seed.getX1();
    float o2 = seed.getX2();
    float o3 = seed.getX3();
    _o1 = o1;
    _o2 = o2;
    _o3 = o3;
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    float[] u = seed.getFaultNormal();
    float[] v = seed.getFaultDipVector();
    float[] w = seed.getFaultStrikeVector();
    updateTransformMap(nu,nv,nw,u,v,w);
    FaultCell[][] cells = new FaultCell[nw][nv];
    // Cell comparator for high-to-low ordering based on fault likelihoods.
    Comparator<FaultCell> flComparator = new Comparator<FaultCell>() {
      public int compare(FaultCell c1, FaultCell c2) {
        if (c1.fl<c2.fl)
          return 1;
        else if (c1.fl>c2.fl)
          return -1;
        else
          return 0;
      }
    };
    float fp = seed.getFp();
    float ft = seed.getFt();
    float fc = seed.getFl();
    seed = new FaultCell(ru,rv,rw,fc,fp,ft);
    cells[rw][rv] = seed;
    // Make a new empty skin.
    FaultSkin skin = new FaultSkin();
    // Make a priority queue of cells, initially with only the seed.
    PriorityQueue<FaultCell> growQueue = 
      new PriorityQueue<FaultCell>(1024,flComparator);
    growQueue.add(seed);
    // While the grow queue is not empty, ...
    while (!growQueue.isEmpty()) {
      FaultCell cell = growQueue.poll();
      if(cell.skin!=null) continue;
      skin.add(cell);
      int uc = cell.getI1();
      int vc = cell.getI2();
      int wc = cell.getI3();
      float vc1 = _vs[0][vc]+o1;
      float vc2 = _vs[1][vc]+o2;
      float vc3 = _vs[2][vc]+o3;
      float wc1 = _ws[0][wc]+o1;
      float wc2 = _ws[1][wc]+o2;
      float wc3 = _ws[2][wc]+o3;
      float xc1 = _us[0][uc]+vc1+wc1-o1;
      float xc2 = _us[1][uc]+vc2+wc2-o2;
      float xc3 = _us[2][uc]+vc3+wc3-o3;
      if(uc<=1||uc>=nu-2) continue;
      if(vc<=1||vc>=nv-2) continue;
      if(wc<=1||wc>=nw-2) continue;
      if(xc1<=1||xc1>=n1-2) continue;
      if(xc2<=1||xc2>=n2-2) continue;
      if(xc3<=1||xc3>=n3-2) continue;
      int c1 = round(xc1);
      int c2 = round(xc2);
      int c3 = round(xc3);
      float fpc = _fp[c3][c2][c1];
      if(_fcgx.findCellsInBox(xc1,xc2,xc3,fpc,2,2,2)) continue;
      int ub = max(uc-5,0);
      int ue = min(uc+5,nu-1);
      int mu = ue-ub+1;
      int dw = _dw;
      float du = 5f;
      if(cell.ca==null) {
        int va = vc-1;
        FaultCell ca = cells[wc][va];
        if(ca!=null) {
          if(abs(ca.x1-cell.x1)<=du)linkAboveBelow(ca,cell);
        } else {
          int ev = max(vc-dw,0);
          int mv = vc-ev+1;
          float[][] fx = new float[mv][mu];
          getArrayAB(-1,ub,ue,vc,wc1,wc2,wc3,fv,fx);
          float[] pik = pick(uc-ub,fx);
          linkAboveBelow(fmin,ub,du,pik,fx,cell,cells,growQueue);
        }
      }
      if(cell.cb==null) {
        int vb = vc+1;
        FaultCell cb = cells[wc][vb];
        if(cb!=null) {
          if(abs(cb.x1-cell.x1)<=du) linkAboveBelow(cell,cb);
        } else {
          int ev = min(vc+dw,nv-1);
          int mv = -vc+ev+1;
          float[][] fx = new float[mv][mu];
          getArrayAB( 1,ub,ue,vc,wc1,wc2,wc3,fv,fx);
          float[] pik = pick(uc-ub,fx);
          linkBelowAbove(fmin,ub,du,pik,fx,cell,cells,growQueue);
        }
      }
      if(cell.cl==null) {
        int wl = wc-1;
        FaultCell cl = cells[wl][vc];
        if(cl!=null) {
         if(abs(cl.x1-cell.x1)<=du) linkLeftRight(cl,cell);
        } else {
          int ew = max(wc-dw,0);
          int mw = wc-ew+1;
          float[][] fx = new float[mw][mu];
          getArrayLR(-1,ub,ue,wc,vc1,vc2,vc3,fv,fx);
          float[] pik = pick(uc-ub,fx);
          linkLeftRight(fmin,ub,du,pik,fx,cell,cells,growQueue);
        }
      }
      if(cell.cr==null) {
        int wr = wc+1;
        FaultCell cr = cells[wr][vc];
        if(cr!=null) {
          if(abs(cr.x1-cell.x1)<=du) linkLeftRight(cell,cr);
        } else {
          int ew = min(wc+dw,nw-1);
          int mw = -wc+ew+1;
          float[][] fx = new float[mw][mu];
          getArrayLR( 1,ub,ue,wc,vc1,vc2,vc3,fv,fx);
          float[] pik = pick(uc-ub,fx);
          linkRightLeft(fmin,ub,du,pik,fx,cell,cells,growQueue);
        }
      }
    }
    FaultSkin skinn = reskin(nw,nw,u,v,w,seed,skin);
    int id = 0;
    for (FaultCell cell:skinn) {
      FaultCell ca = cell.ca;
      FaultCell cb = cell.cb;
      FaultCell cl = cell.cl;
      FaultCell cr = cell.cr;
      if(ca!=null&&ca.skin==null) cell.ca=null;
      if(cb!=null&&cb.skin==null) cell.cb=null;
      if(cl!=null&&cl.skin==null) cell.cl=null;
      if(cr!=null&&cr.skin==null) cell.cr=null;
      backToXyz(ru,u,cell,fv);
      cell.id = id;
      id++;
    }
    return skinn;
  }


  public FaultSkin reskin(
    int nv, int nw, float[] u, float[] v, float[] w,
    FaultCell seed, FaultSkin skin) {
    float[][] ux = new float[nw][nv];
    float[][] fx = new float[nw][nv];
    int ve = 0;
    int we = 0;
    int vb = nv;
    int wb = nw;
    for (FaultCell cell:skin) {
      int iv = cell.i2;
      int iw = cell.i3;
      ux[iw][iv] = cell.x1;
      fx[iw][iv] = cell.fl;
      if(iv<vb) vb=iv;
      if(iw<wb) wb=iw;
      if(iv>ve) ve=iv;
      if(iw>we) we=iw;
    }
    int mv = ve-vb+1;
    int mw = we-wb+1;
    float[][] fxs = copy(mv,mw,vb,wb,fx);
    float[][] uxs = copy(mv,mw,vb,wb,ux);
    float[][] uss = smooth(4,4,fxs,uxs);
    float[][] fss = new float[mw][mv];
    _rgf.apply00(fxs,fss);
    float du = 5f;
    float fmin = 0.2f;
    int vs = seed.i2;
    int ws = seed.i3;
    float xvs = seed.x2;
    float xws = seed.x3;
    float xus = uss[ws-wb][vs-vb];
    float fls = fss[ws-wb][vs-vb];
    float[][][] fpt = surfaceStrikeAndDip(u,v,w,uss);
    float fps = fpt[0][ws-wb][vs-vb];
    float fts = fpt[1][ws-wb][vs-vb];
    FaultCell[][] cells = new FaultCell[mw][mv];
    seed = new FaultCell(xus,xvs,xws,fls,fps,fts);
    FaultSkin skinn = new FaultSkin();
    // Cell comparator for high-to-low ordering based on fault likelihoods.
    Comparator<FaultCell> flComparator = new Comparator<FaultCell>() {
      public int compare(FaultCell c1, FaultCell c2) {
        if (c1.fl<c2.fl)
          return 1;
        else if (c1.fl>c2.fl)
          return -1;
        else
          return 0;
      }
    };
    // Make a priority queue of cells, initially with only the seed.
    PriorityQueue<FaultCell> growQueue = 
      new PriorityQueue<FaultCell>(1024,flComparator);
    growQueue.add(seed);
    // While the grow queue is not empty, ...
    while (!growQueue.isEmpty()) {
      FaultCell cell = growQueue.poll();
      if(cell.skin!=null) continue;
      skinn.add(cell);
      int uc = cell.i1;
      int vc = cell.i2;
      int wc = cell.i3;
      int vm = vc-vb;
      int wm = wc-wb;
      if(vm<=1||vm>=mv-2) continue;
      if(wm<=1||wm>=mw-2) continue;
      if(cell.ca==null) {
        FaultCell ca = cells[wm][vm-1];
        if(ca != null) {
          if(abs(ca.x1-uc)<du) linkAboveBelow(ca,cell);
        } else {
          float fla = fss[wm][vm-1];
          float xau = uss[wm][vm-1];
          float fpa = fpt[0][wm][vm-1];
          float fta = fpt[1][wm][vm-1];
          if(fla>fmin&&abs(uc-xau)<du) {
            ca = new FaultCell(xau,vc-1,wc,fla,fpa,fta);
            linkAboveBelow(ca,cell);
            cells[wm][vm-1] = ca;
            growQueue.add(ca);
          }
        }
      }

      if(cell.cb==null) {
        FaultCell cb = cells[wm][vm+1];
        if(cb != null) {
          if(abs(cb.x1-uc)<du) linkAboveBelow(cell,cb);
        } else {
          float flb = fss[wm][vm+1];
          float xbu = uss[wm][vm+1];
          float fpb = fpt[0][wm][vm+1];
          float ftb = fpt[1][wm][vm+1];
          if(flb>fmin&&abs(uc-xbu)<du) {
            cb = new FaultCell(xbu,vc+1,wc,flb,fpb,ftb);
            linkAboveBelow(cell,cb);
            cells[wm][vm+1] = cb;
            growQueue.add(cb);
          }
        }
      }

      if(cell.cl==null) {
        FaultCell cl = cells[wm-1][vm];
        if(cl != null) {
          if(abs(cl.x1-uc)<du) linkLeftRight(cl,cell);
        } else {
        float fll = fss[wm-1][vm];
        float xlu = uss[wm-1][vm];
        float fpl = fpt[0][wm-1][vm];
        float ftl = fpt[1][wm-1][vm];
        if(fll>fmin&&abs(uc-xlu)<du) {
          cl = new FaultCell(xlu,vc,wc-1,fll,fpl,ftl);
          linkLeftRight(cl,cell);
          growQueue.add(cl);
          cells[wm-1][vm] = cl;
        }
        }
      }

      if(cell.cr==null) {
        FaultCell cr = cells[wm+1][vm];
        if(cr != null) {
          if(abs(cr.x1-uc)<du) linkLeftRight(cell,cr);
        } else {
        float flr = fss[wm+1][vm];
        float xru = uss[wm+1][vm];
        float fpr = fpt[0][wm+1][vm];
        float ftr = fpt[1][wm+1][vm];
        if(flr>fmin&&abs(uc-xru)<du) {
          cr = new FaultCell(xru,vc,wc+1,flr,fpr,ftr);
          linkLeftRight(cell,cr);
          growQueue.add(cr);
          cells[wm+1][vm] = cr;
        }
        }
      }
    }
    return skinn;
  }

  private float[][][] surfaceStrikeAndDip(
    float[] u, float[] v, float[] w, float[][] sfu) {
    int nw = sfu.length;
    int nv = sfu[0].length;
    float[][] sf1 = new float[nw][nv];
    float[][] sf2 = new float[nw][nv];
    float[][] fp  = new float[nw][nv];
    float[][] ft  = new float[nw][nv];
    for (int iw=0; iw<nw; ++iw) {
    for (int iv=0; iv<nv; ++iv) {
      int mw = iw-1; mw = max(mw,0);
      int mv = iv-1; mv = max(mv,0);
      int pw = iw+1; pw = min(pw,nw-1);
      int pv = iv+1; pv = min(pv,nv-1);
      int dw = pw-mw;
      int dv = pv-mv;
      sf1[iw][iv] = (sfu[iw][pv]-sfu[iw][mv])/dv;
      sf2[iw][iv] = (sfu[pw][iv]-sfu[mw][iv])/dw;
      float a1 = -1f;
      float a2 = sf1[iw][iv];
      float a3 = sf2[iw][iv];
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
      ft[iw][iv] = toDegrees(acos(-t1));
      fp[iw][iv] = range360(toDegrees(atan2(-t3,t2)));
    }}
    return new float[][][]{fp,ft};
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



  public float[][] smooth(
    float sig1, float sig2, float[][] w, float[][] x) {
    int n3 = w.length;
    int n2 = w[0].length;
    float[][] b = new float[n3][n2];
    float[][] r = new float[n3][n2];
    makeRhs(x,w,b);
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(sig1,sig2);
    A2 a2 = new A2(smoother2,w);
    CgSolver cs = new CgSolver(0.001,200);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(r);
    return vr.getArray();
  }

  private void makeRhs(
    float[][] x, float[][] w, float[][] b) {
    int n3 = x.length;
    int n2 = x[0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      float wi = w[i3][i2];
      b[i3][i2] = x[i3][i2]*wi*wi;
    }}
  }

    // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] wp) 
    {
      _s2 = s2;
      _wp = wp;
      float n2 = wp.length;
      float n1 = wp[0].length;
      _sc = 4f*sum(wp)/(n1*n2);
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      v2y.zero();
      _s2.apply(z);
      addAndScale(-_sc,z,y);
      applyLhs(_wp,z,y);
      _s2.applyTranspose(y);
      addAndScale( _sc,x,y);
    }
    private float _sc;
    private Smoother2 _s2;
    private float[][] _wp;
  }

  private static void applyLhs(float[] wp, float[] x, float[] y) {
    int n1 = wp.length;
    for (int i1=0; i1<n1; ++i1)
      y[i1] += wp[i1]*wp[i1]*x[i1];
  }


  private static void applyLhs(float[][] wp, float[][] x, float[][] y) {
    int n2 = wp.length;
    for (int i2=0; i2<n2; ++i2)
      applyLhs(wp[i2],x[i2],y[i2]);
  }

  private static void addAndScale(float sc, float[][] x, float[][] y) {
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      addAndScale(sc,x[i2],y[i2]);
    }
  }

  private static void addAndScale(float sc, float[] x, float[] y) {
    int n1 = x.length;
    for (int i1=0; i1<n1; ++i1) {
      y[i1] += sc*x[i1];
    }
  }


  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(float sigma1, float sigma2) {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
    }
    public void apply(float[][] x) {
      smooth2(_sigma2,x);
      smooth1(_sigma1,x);
    }
    public void applyTranspose(float[][] x) {
      smooth1(_sigma1,x);
      smooth2(_sigma2,x);
    }
    private float _sigma1,_sigma2;
  }


  // Smoothing for dimension 2.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] x) {
    if (sigma<1.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply2(x,x);
  }

  private void backToXyz(
    int ru, float[] u, FaultCell cell, float[][][] fv) {
    int iv = cell.getI2();
    int iw = cell.getI3();
    float xu = cell.getX1()-ru;
    float dw1 = _ws[0][iw]+_o1;
    float dw2 = _ws[1][iw]+_o2;
    float dw3 = _ws[2][iw]+_o3;
    float x1 = (u[0]*xu+_vs[0][iv]+dw1);
    float x2 = (u[1]*xu+_vs[1][iv]+dw2);
    float x3 = (u[2]*xu+_vs[2][iv]+dw3);
    x1 = max(x1,0);
    x2 = max(x2,0);
    x3 = max(x3,0);
    x1 = min(x1,_n1-1);
    x2 = min(x2,_n2-1);
    x3 = min(x3,_n3-1);
    cell.x1 = x1;
    cell.x2 = x2;
    cell.x3 = x3;
    int c1 = round(x1);
    int c2 = round(x2);
    int c3 = round(x3);
    cell.i1 = c1;
    cell.i2 = c2;
    cell.i3 = c3;
    cell.fl = fv[c3][c2][c1];
  }

  private void linkAboveBelow(
    float fmin, float ub, float du,
    float[] pik, float[][] fx, FaultCell cell, 
    FaultCell[][] cells, PriorityQueue<FaultCell> growQueue) {
    int vc = cell.i2;
    int wc = cell.i3;
    float uf = cell.x1;
    float ft = cell.ft;
    int mv = pik.length;
    FaultCell cb = cell;
    for (int iv=1; iv<mv; iv++) {
      float up = pik[iv];
      int iu = round(up);
      float fi = fx[iv][iu];
      int va = vc-iv;
      float pi = getFp(up+ub,va,wc);
      float dp = abs(cb.fp-pi);
      if(dp>180) dp = 360-dp;
      if(fi<fmin||dp>_dpmax) break;
      FaultCell ct = cells[wc][va];
      if(ct==null) {
        up += ub;
        FaultCell ca = new FaultCell(up,va,wc,fi,pi,ft);
        growQueue.add(ca);
        cells[wc][va] = ca;
        linkAboveBelow(ca,cb);
        cb = ca;
      } else if(ct.cb==null&&abs(ct.x1-uf)<=du) {
        linkAboveBelow(ct,cb);
        cb = ct;
      } else { break;}
    }
  }

  private float getFp(float u, float v, float w) {
    int iu = round(u);
    int iv = round(v);
    int iw = round(w);
    int i1 = round(_us[0][iu]+_vs[0][iv]+_ws[0][iw]+_o1);
    int i2 = round(_us[1][iu]+_vs[1][iv]+_ws[1][iw]+_o2);
    int i3 = round(_us[2][iu]+_vs[2][iv]+_ws[2][iw]+_o3);
    i1 = max(i1,0);
    i2 = max(i2,0);
    i3 = max(i3,0);
    i1 = min(i1,_n1-1);
    i2 = min(i2,_n2-1);
    i3 = min(i3,_n3-1);
    return _fp[i3][i2][i1];
  }


  private void linkBelowAbove(
    float fmin, float ub, float du,
    float[] pik, float[][] fx, FaultCell cell, 
    FaultCell[][] cells, PriorityQueue<FaultCell> growQueue) {
    int vc = cell.i2;
    int wc = cell.i3;
    float uf = cell.x1;
    float ft = cell.ft;
    int mv = pik.length;
    FaultCell ca = cell;
    for (int iv=1; iv<mv; iv++) {
      float up = pik[iv];
      int iu = round(up);
      float fi = fx[iv][iu];
      int vb = vc+iv;
      float pi = getFp(up+ub,vb,wc);
      float dp = abs(ca.fp-pi);
      if(dp>180) dp = 360-dp;
      if(fi<fmin||dp>_dpmax) break;
      FaultCell ct = cells[wc][vb];
      if(ct==null) {
        up += ub;
        FaultCell cb = new FaultCell(up,vb,wc,fi,pi,ft);
        growQueue.add(cb);
        cells[wc][vb] = cb;
        linkAboveBelow(ca,cb);
        ca = cb;
      } else if(ct.ca==null&&abs(ct.x1-uf)<=du) {
        linkAboveBelow(ca,ct);
        ca = ct;
      } else { break;}
    }
  }

  private void linkLeftRight(
    float fmin, float ub, float du,
    float[] pik, float[][] fx, FaultCell cell, 
    FaultCell[][] cells, PriorityQueue<FaultCell> growQueue) {
    int vc = cell.i2;
    int wc = cell.i3;
    float uf = cell.x1;
    float ft = cell.ft;
    int mw = pik.length;
    FaultCell cr = cell;
    for (int iw=1; iw<mw; iw++) {
      float up = pik[iw];
      int iu = round(up);
      float fi = fx[iw][iu];
      int wl = wc-iw;
      float pi = getFp(up+ub,vc,wl);
      float dp = abs(cr.fp-pi);
      if(dp>180) dp = 360-dp;
      if(fi<fmin||dp>_dpmax) break;
      FaultCell ct = cells[wl][vc];
      if(ct==null) {
        up += ub;
        FaultCell cl = new FaultCell(up,vc,wl,fi,pi,ft);
        growQueue.add(cl);
        cells[wl][vc] = cl;
        linkLeftRight(cl,cr);
        cr = cl;
      } else if(ct.cr==null&&abs(ct.x1-uf)<=du) {
        linkLeftRight(ct,cr);
        cr = ct;
      } else { break;}
    }
  }

  private void linkRightLeft(
    float fmin, float ub, float du,
    float[] pik, float[][] fx, FaultCell cell, 
    FaultCell[][] cells, PriorityQueue<FaultCell> growQueue) {
    int vc = cell.i2;
    int wc = cell.i3;
    float uf = cell.x1;
    float ft = cell.ft;
    int mw = pik.length;
    FaultCell cl = cell;
    for (int iw=1; iw<mw; iw++) {
      float up = pik[iw];
      int iu = round(up);
      float fi = fx[iw][iu];
      int wr = wc+iw;
      float pi = getFp(up+ub,vc,wr);
      float dp = abs(cl.fp-pi);
      if(dp>180) dp = 360-dp;
      if(fi<fmin||dp>_dpmax) break;
      FaultCell ct = cells[wr][vc];
      if(ct==null) {
        up += ub;
        FaultCell cr = new FaultCell(up,vc,wr,fi,pi,ft);
        growQueue.add(cr);
        cells[wr][vc] = cr;
        linkLeftRight(cl,cr);
        cl = cr;
      } else if(ct.cl==null&&abs(ct.x1-uf)<=du) {
        linkLeftRight(cl,ct);
        cl = ct;
      } else { break;}
    }
  }

  private void getArrayLR(
    int dir, int ub, int ue,
    int wc, float vc1, float vc2, float vc3,
    float[][][] fv, float[][] fx) {
    int mw = fx.length;
    int n3 = fv.length;
    int n2 = fv[0].length;
    int n1 = fv[0][0].length;
    for (int iw=0; iw<mw; iw++) {
      int wt = wc;
      if(dir==-1) wt -= iw;
      else        wt += iw;
      float wx1 = _ws[0][wt];
      float wx2 = _ws[1][wt];
      float wx3 = _ws[2][wt];
      for (int iu=ub; iu<=ue; iu++) { 
        int i1 = round(_us[0][iu]+vc1+wx1);
        int i2 = round(_us[1][iu]+vc2+wx2);
        int i3 = round(_us[2][iu]+vc3+wx3);
        if(i1<0||i1>=n1) continue;
        if(i2<0||i2>=n2) continue;
        if(i3<0||i3>=n3) continue;
        fx[iw][iu-ub] = fv[i3][i2][i1];
      }
    }
  }

  private void getArrayAB(
    int dir, int ub, int ue,
    int vc, float wc1, float wc2, float wc3,
    float[][][] fv, float[][] fx) {
    int mv = fx.length;
    int n3 = fv.length;
    int n2 = fv[0].length;
    int n1 = fv[0][0].length;
    for (int iv=0; iv<mv; iv++) {
      int vt = vc;
      if(dir==-1) vt -= iv;
      else        vt += iv;
      float vx1 = _vs[0][vt];
      float vx2 = _vs[1][vt];
      float vx3 = _vs[2][vt];
    for (int iu=ub; iu<=ue; iu++) { 
      int i1 = round(_us[0][iu]+vx1+wc1);
      int i2 = round(_us[1][iu]+vx2+wc2);
      int i3 = round(_us[2][iu]+vx3+wc3);
      if(i1<0||i1>=n1) continue;
      if(i2<0||i2>=n2) continue;
      if(i3<0||i3>=n3) continue;
      fx[iv][iu-ub] = fv[i3][i2][i1];
    }}
  }



  private float[] pick(int i0, float[][] fx) {
    OptimalPathPicker opp = new OptimalPathPicker(4,_an);
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = opp.applyTransform(fx);
    float[][] wh = opp.applyForWeight(ft);
    float[][] tm = zerofloat(n2,n1);
    float[] pk = opp.forwardPick(i0,wh,tm);
    return pk;
  }


  // Methods to link mutually nabors.
  private void linkAboveBelow(FaultCell ca, FaultCell cb) {
    cb.ca = ca;
    ca.cb = cb;
  }
  private void linkLeftRight(FaultCell cl, FaultCell cr) {
    cr.cl = cl;
    cl.cr = cr;
  }

  private void updateTransformMap(
    int nu, int nv, int nw,
    float[] u, float[] v, float[] w) {
    _us = new float[3][nu];
    _vs = new float[3][nv];
    _ws = new float[3][nw];
    updateTransformMap(u,_us);
    updateTransformMap(v,_vs);
    updateTransformMap(w,_ws);
  }

  private void updateTransformMap(
    float[] u, float[][] us) {
    int nr = us[0].length;
    int r = (nr-1)/2;
    for (int i=1; i<=r; ++i) {
      int kp = r+i;
      int km = r-i;
      float iu1 = i*u[0];
      float iu2 = i*u[1];
      float iu3 = i*u[2];
      us[0][kp] =  iu1;
      us[1][kp] =  iu2;
      us[2][kp] =  iu3;
      us[0][km] = -iu1;
      us[1][km] = -iu2;
      us[2][km] = -iu3;
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  // private
  private float[][] _us,_vs,_ws;
  private float[][][] _fp, _ft;
  private float _o1,_o2,_o3;
  private int _n1,_n2,_n3;
  private FaultCellGrid _fcgx;
  private float _dpmax = 180;
  private int _dw=10;
  private float _an=0.2f;
  private SincInterpolator _si = new SincInterpolator();
  private RecursiveGaussianFilterP _rgf = new RecursiveGaussianFilterP(8);
  private RecursiveGaussianFilterP _rgfG = new RecursiveGaussianFilterP(1);
  private RecursiveGaussianFilterP _rgfu = new RecursiveGaussianFilterP(2);
  private RecursiveGaussianFilterP _rgfv = new RecursiveGaussianFilterP(4);
  private RecursiveGaussianFilterP _rgfw = new RecursiveGaussianFilterP(4);
}
