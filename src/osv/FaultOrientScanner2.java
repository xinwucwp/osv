/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package osv;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.*;

/**
 * Scan for fault orientations.
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.07.15
 */

public class FaultOrientScanner2 {
  /**
   * Constructs a scanner with specified parameters.
   * @param sigmaTheta half-width for smoothing along dip of edges.
   */
  public FaultOrientScanner2(double sigma1) {
    _sigma1 = (float)sigma1;
    _si = new SincInterpolator();
    _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
  }

  /**
   * Gets a sampling of edge dip theta appropriate for this scanner.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   */
  public Sampling getThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma1,thetaMin,thetaMax);
  }

  /**
   * Scans a specified image for edge dips.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   * @param g the image to be scanned.
   * @return array {el,et} of edge gradients and dips.
   */
  public float[][][] scan(double thetaMin, double thetaMax, float[][] g) {
    Sampling st = makeThetaSampling(thetaMin,thetaMax);
    //Sampling st = new Sampling(18,10,0);
    return scanTheta(st,g);
  }

  /**
   * Scans a specified image for edge dips.
   * @param thetaMin minimum edge dip, in degrees.
   * @param thetaMax maximum edge dip, in degrees.
   * @param g the image to be scanned.
   * @return array {el,et} of edge gradients and dips.
   */
  public float[][][] scanDip(double thetaMin, double thetaMax, float[][] g) {
    int n2 = g.length;
    int n1 = g[0].length;
    Sampling st1 = makeThetaSampling(90-thetaMax,90-thetaMin);
    Sampling st2 = makeThetaSampling(90+thetaMin,90+thetaMax);
    float[][][] fp1 = scanTheta(st1,g);
    float[][][] fp2 = scanTheta(st2,g);
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fx1 = fp1[0][i2][i1];
      float fx2 = fp2[0][i2][i1];
      if(fx1<fx2) {
        fp1[0][i2][i1] = fx2;
        fp1[1][i2][i1] = fp2[1][i2][i1];
      }
    }}
    return fp1;
  }


  public float[][] striveConvert(float[][] fp) {
    int n2 = fp.length;
    int n1 = fp[0].length;
    float[][] fc = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float fpi = 90-fp[i2][i1];
      if(fpi<0f) fpi+=180;
      fc[i2][i1] = fpi;
    }}
    return fc;
  }

  /**
   * Thins fault images to include only ridges in image edges.
   * After thinning, may be only one voxel wide. Thinned fault strikes and
   * dips are set to zero where thinned fault likelihoods are zero.
   * @param fet array {el,et} of edge gradeitns and dips.
   * @return array {elt,ett} of thinned edge gradeints and dips.
   */
  public float[][][] thin(float[][][] fet) {
    int n2 = fet[0].length;
    final int n1 = fet[0][0].length;
    final float[][] f = fet[0];
    final float[][] t = fet[1];
    final float[][] ff = new float[n2][n1];
    final float[][] tt = new float[n2][n1];
    final float pi = (float)(Math.PI/180.0);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float ti = t[i2][i1]*pi;
        float d1 = sin(ti);
        float d2 = cos(ti);
        float x1p = i1+d1;
        float x2p = i2+d2;
        float x1m = i1-d1;
        float x2m = i2-d2;
        float fi = f[i2][i1];
        float fp = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1p,x2p);
        float fm = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,f,x1m,x2m);
        if(fp<fi&&fm<fi) {
          ff[i2][i1] = fi;
          tt[i2][i1] = t[i2][i1];
        }
      }
    }});
    return new float[][][]{ff,tt};
  }

  public float[][] edgeLikeFit2(int r, float[][] fl) {
    int n2 = fl.length;
    int n1 = fl[0].length;
    float[][] flr = new float[n2][n1];
    for (int i2=r; i2<n2-r-1; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] ft = new float[r*2+1];
      for (int r2=-r;r2<=r;r2++)
        ft[r2+r] = fl[i2+r2][i1];
      float[] abc = parabolicFit(ft);
      float a = abc[0];
      float b = abc[1];
      float s = -abs(2*a*r+b)/(2*a);
      flr[i2][i1] = fl[i2][i1]/(s+0.0001f);
    }}
    return flr;
  }

  public float[][] edgeLikeFit(int r, float[][] el, float[][] et) {
    int n2 = el.length;
    int n1 = el[0].length;
    float[][] elr = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float[] ft = new float[r*2+1];
      float ti = et[i2][i1]+90f;
      if(ti>180) ti -= 180;
      ti *= pi;
      float t1 = cos(ti);
      float t2 = sin(ti);
      for (int ir=-r;ir<=r;ir++)
        ft[ir+r] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,el,i1+ir*t1,i2+ir*t2);
      float[] abc = parabolicFit(ft);
      float a = abc[0];
      float b = abc[1];
      float s = -abs(2*a*r+b)/(2*a);
      elr[i2][i1] = el[i2][i1]/(s+0.0001f);
    }}
    return elr;
  }


  public float[][] edgeLikeFitG(int r, float[][] el, float[][] et) {
    int n2 = el.length;
    int n1 = el[0].length;
    float[][] elr = new float[n2][n1];
    float pi = (float)(Math.PI/180.0);
    SincInterpolator si = new SincInterpolator();
    RecursiveGaussianFilterP rgf = new RecursiveGaussianFilterP(r);
    float[] ft = new float[10*r+1];
    float[] fs = new float[10*r+1];
    float[] f1 = new float[10*r+1];
    float[] f2 = new float[10*r+1];
    int hr = 5*r;
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      float ti = et[i2][i1]+90f;
      if(ti>180) ti -= 180;
      ti *= pi;
      float t1 = cos(ti);
      float t2 = sin(ti);
      for (int ir=-hr;ir<=hr;ir++)
        ft[ir+hr] = si.interpolate(n1,1.0,0.0,n2,1.0,0.0,el,i1+ir*t1,i2+ir*t2);
      rgf.apply0(ft,fs);
      rgf.apply1(ft,f1);
      rgf.apply2(ft,f2);
      elr[i2][i1] = fs[hr]*(-f2[hr]/(abs(f1[hr])+0.0001f));
    }}
    return elr;
  }

  public float[] parabolicFit(float[] f) {
    int n1 = f.length;
    double[][] A = new double[3][3];
    double[][] B = new double[3][1];
    double s4=0.0;
    double s3=0.0;
    double s2=0.0;
    double s1=0.0;
    double s0=0.0;
    double b0=0.0;
    double b1=0.0;
    double b2=0.0;
    for (int i1=0; i1<n1; ++i1) {
      double x1 = i1;
      double x2 = i1*x1;
      double x3 = i1*x2;
      double x4 = i1*x3;
      s0 += 1.;
      s1 += x1;
      s2 += x2;
      s3 += x3;
      s4 += x4;
      b0 += x2*f[i1];
      b1 += x1*f[i1];
      b2 += f[i1];
    }
    A[0][0] = s4;
    A[1][0] = s3;
    A[2][0] = s2;
    A[0][1] = s3;
    A[1][1] = s2;
    A[2][1] = s1;
    A[0][2] = s2;
    A[1][2] = s1;
    A[2][2] = s0;
    B[0][0] = b0;
    B[1][0] = b1;
    B[2][0] = b2;
    DMatrix da = new DMatrix(A);
    DMatrix db = new DMatrix(B);
    DMatrix dx = da.solve(db);
    double[][] x = dx.get();
    float a = (float)x[0][0];
    float b = (float)x[1][0];
    float c = (float)x[2][0];
    return new float[]{a,b,c};
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1;
  private SincInterpolator _si;

  private static final float NO_DIP    = -0.00001f;

  private Sampling makeThetaSampling(double thetaMin, double thetaMax) {
    return angleSampling(_sigma1,thetaMin,thetaMax);
  }
  private static Sampling angleSampling(
    double sigma, double amin, double amax)
  {
    double fa = amin;
    double da = toDegrees(0.5/sigma);
    int na = 1+(int)((amax-amin)/da);
    da = (amax>amin)?(amax-amin)/(na-1):1.0;
    return new Sampling(na,da,fa);
  }

  // Smoothing filter
  private static RecursiveExponentialFilter makeRef(double sigma) {
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
    return ref;
  }

  // Scans over all edge dips theta. 
  private float[][][] scanTheta(Sampling thetaSampling, float[][] g) {
    final int n2 = g.length;
    final int n1 = g[0].length;
    final float[][] f = new float[n2][n1];
    final float[][] t = new float[n2][n1];
    final SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    int nt = thetaSampling.getCount();
    float[][][] rsc = getRotationAngleMap(n1,n2);
    RecursiveExponentialFilter ref2 = makeRef(1);
    RecursiveExponentialFilter ref1 = makeRef(_sigma1);
    for (int it=0; it<nt; ++it) {
      System.out.println(it+"/"+(nt-1)+" done...");
      float ti = (float)thetaSampling.getValue(it);
      float theta = toRadians(ti);
      float[][] gr = rotate(theta,rsc,g);
      ref2.apply2(gr,gr);
      ref1.apply1(gr,gr);
      float[][] s2 = unrotate(n1,n2,theta,rsc,gr);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float st = s2[i2][i1];
        st *= st;
        st *= st;
        st = 1-st;
        if (st>f[i2][i1]) {
          f[i2][i1] = st;
          t[i2][i1] = ti;
        }
      }}
    }
    return new float[][][]{f,t};
  }

  public float[][] rotate(float theta, 
    float[][][] rsc, float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    final int nr = rsc[0].length;
    final int ri = (nr-1)/2;
    final int h2 = (int)floor(n2/2);
    final int h1 = (int)floor(n1/2);
    int r21 = abs(round(h2*cos(theta)+h1*sin(theta)));
    int r22 = abs(round(h2*cos(theta)-h1*sin(theta)));
    int r11 = abs(round(h1*cos(theta)-h2*sin(theta)));
    int r12 = abs(round(h1*cos(theta)+h2*sin(theta)));
    int r2 = max(r21,r22);
    int r1 = max(r11,r12);
    final int m2 = r2*2+1;
    final int m1 = r1*2+1;
    final float st = sin(theta);
    final float ct = cos(theta);
    final float[][] fr = new float[m2][m1];
    Parallel.loop(m2,new Parallel.LoopInt() {
    public void compute(int i2) {
      int k2 = i2+ri-r2;
      float[] rr2 = rsc[0][k2];
      float[] ss2 = rsc[1][k2];
      float[] cc2 = rsc[2][k2];
      for (int i1=0; i1<m1; ++i1) {
        int k1 = i1+ri-r1;
        float rk = rr2[k1]; 
        float sk = ss2[k1]; 
        float ck = cc2[k1]; 
        float sp = sk*ct-ck*st;
        float cp = ck*ct+sk*st;
        float x1 = rk*cp+h1;
        float x2 = rk*sp+h2;
        fr[i2][i1] = _si.interpolate(n1,1.0,0.0,n2,1.0,0.0,fx,x1,x2);
      }
    }});
    return fr;
  }

  // unrotate a 3D iamge around axis1
  public float[][] unrotate(int n1, int n2, float theta, 
    float[][][] rsc, float[][] fr) {
    final int h2 = (int)floor(n2/2);
    final int h1 = (int)floor(n1/2);
    final int nr = rsc[0].length;
    final int ri = (nr-1)/2;
    final int m2 = fr.length;
    final int m1 = fr[0].length;
    final int r1 = (m1-1)/2;
    final int r2 = (m2-1)/2;
    final float st = sin(theta);
    final float ct = cos(theta);
    final float[][] fu = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      int k2 = i2+ri-h2;
      float[] rr2 = rsc[0][k2];
      float[] ss2 = rsc[1][k2];
      float[] cc2 = rsc[2][k2];
      for (int i1=0; i1<n1; ++i1) {
        int k1 = i1+ri-h1;
        float rk = rr2[k1]; 
        float sk = ss2[k1]; 
        float ck = cc2[k1]; 
        float sp = sk*ct+ck*st;
        float cp = ck*ct-sk*st;
        float x1 = rk*cp+r1;
        float x2 = rk*sp+r2;
        fu[i2][i1] = _si.interpolate(m1,1.0,0.0,m2,1.0,0.0,fr,x1,x2);
      }
    }});
    return fu;
  }

  public float[][][] getRotationAngleMap(int n1, int n2) {
    int h1 = (int)ceil(n1/2);
    int h2 = (int)ceil(n2/2);
    int hr = round(sqrt(h1*h1+h2*h2));
    int nr = 2*hr+1;
    float[][] sr = new float[nr][nr];
    float[][] cr = new float[nr][nr];
    float[][] rr = new float[nr][nr];
    for (int k2=-hr; k2<=hr; ++k2) {
    for (int k1=-hr; k1<=hr; ++k1) {
      int i1 = k1+hr;
      int i2 = k2+hr;
      float ri = sqrt(k1*k1+k2*k2);
      if(ri==0f) continue;
      float rc = 1.0f/ri;
      rr[i2][i1] = ri;
      sr[i2][i1] = k2*rc;
      cr[i2][i1] = k1*rc;
    }}
    return new float[][][]{rr,sr,cr};
  }

}
