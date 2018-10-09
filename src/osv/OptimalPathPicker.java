package osv;

import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Optimal path picking
 * @author Xinming Wu 
 * @version 2016.03.23
 */

public class OptimalPathPicker {
 
  public OptimalPathPicker(int gate, float an) {
    _gate = gate;
    _an = an;
  }


  public float[][] applyTransform(float[][] fx) {
    int n2 = fx.length;
    int n1 = fx[0].length;
    float[][] ft = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      ft[i1][i2] = fx[i2][i1];
    }}
    return ft;
  }
  public float[][] applyForWeight(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i1][i2]  = exp(-vel[i2][i1]);
    }}
    return w;
  }

  public float[][] applyForWeightX(float[][] vel) {
    int n2 = vel.length;
    int n1 = vel[0].length;
    float[][] w = new float[n1][n2];
    for (int i2=0; i2<n2; ++i2) {
    for (int i1=0; i1<n1; ++i1) {
      w[i2][i1]  = exp(-vel[i2][i1]);
    }}
    return w;
  }


  public float[][] accumulateInline(final float[][][] vel) {
    final int n3 = vel.length;
    final int n2 = vel[0].length;
    final int n1 = vel[0][0].length;
    final float[][] p1 = new float[n2][n3];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] w1 = new float[n3][n1];
      float[][] tf = new float[n3][n1];
      float[][] tb = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(n1-1,w1,tf);
      int i0 = round(p1[i2][n3-1]);
      i0 = min(i0,n1-1); i0 = max(i0,0);
      p1[i2] = backwardPick(i0,w1,tb);
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        vel[i3][i2][i1] = tf[i3][i1]+tb[i3][i1];
      }}
    }});
    return p1;
  }

  public float[][] accumulateInline(final int i10, final float[][][] vel) {
    final int n3 = vel.length;
    final int n2 = vel[0].length;
    final int n1 = vel[0][0].length;
    final float[][] p1 = new float[n2][n3];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] w1 = new float[n3][n1];
      float[][] tf = new float[n3][n1];
      float[][] tb = new float[n3][n1];
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        w1[i3][i1] = exp(-vel[i3][i2][i1]);
      }}
      p1[i2] = forwardPick(i10,w1,tf);
      int i0 = round(p1[i2][n3-1]);
      i0 = min(i0,n1-1); i0 = max(i0,0);
      p1[i2] = backwardPick(i10,w1,tb);
      for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        vel[i3][i2][i1] = tf[i3][i1]+tb[i3][i1];
      }}
    }});
    return p1;
  }


  public float[][] accumulateCrossline(final float[] p, final float[][][] vel){
    final int n3 = vel.length;
    final int n2 = vel[0].length;
    final int n1 = vel[0][0].length;
    final float[][] p2 = new float[n3][n2];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] vel3 = vel[i3];
      float[][] tf = new float[n2][n1];
      float[][] tb = new float[n2][n1];
      float[][] w2 = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        w2[i2][i1] = exp(-vel3[i2][i1]);
      }}
      int i0 = round(p[i3]);
      i0 = min(i0,n1-1);i0 = max(i0,0);
      p2[i3] = forwardPick(i0,w2,tf);
      i0 = round(p2[i3][n2-1]);
      i0 = min(i0,n1-1); i0 = max(i0,0);
      p2[i3] = backwardPick(i0,w2,tb);
      for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        vel3[i2][i1] = tf[i2][i1]+tb[i2][i1];
      }}
    }});
    return p2;
  }


  public float[] backwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    float[] prob = new float[_gate*2-1];
    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    return p;
  }



  public float[] backwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-1][i1]+wx[n2-1][i0]);
	    tt[n2-1][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[n2-2][i1]+wx[n2-1][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[n2-2][i1] = i0;
	    tt[n2-2][i1]=prev[i1];
    }

    float[] prob = new float[_gate*2-1];
    for (int i2=n2-3; i2>=0; i2--) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2+1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    forwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      //tx[i2][i1]  = tt[i2][i1];
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }

  public float[] forwardPick(int i0, float[][] wx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    float[] prob = new float[_gate*2-1];
    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    return p;
  }


  public float[] forwardPick(int i0, float[][] wx, float[][] tx) {
    int n2 = wx.length;
    int n1 = wx[0].length;
    float[] p = new float[n2];
    float[][] tt = new float[n2][n1];
    float[][] what = new float[n2][n1];
    float[] prev = new float[n1];
    float[] next = new float[n1];
    float[] dist = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      dist[i1] = sqrt(i1*i1+_an*_an); 
    }
	  for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[0][i1]+wx[0][i0]);
	    tt[0][i1] = abs(i1-i0)*wi;
	  }
    
    for (int i1=0; i1<n1; i1++) {
	    float wi = 0.5f*(wx[1][i1]+wx[0][i0]);
	    prev[i1] = dist[abs(i1-i0)]*wi;
	    what[1][i1] = i0;
	    tt[1][i1]=prev[i1];
    }

    float[] prob = new float[_gate*2-1];
    for (int i2=2; i2<n2; i2++) {
	    for (int i1=0; i1<n1; i1++) {
	      float wi = wx[i2][i1];
	      int ib = max(i1-_gate,-1);
	      int ie = min(i1+_gate,n1);
	      float c = FLT_MAX;
	      int ic = -1;
	      for (int i=ib+1; i<ie; i++) {
		      float w2 = 0.5f*(wi+wx[i2-1][i]);
		      float d = dist[abs(i1-i)]*w2+prev[i];
		      int it = i-ib-1;
		      if (d < c) {
		        c =	d;
		        ic = it;
		      }
		      prob[it]=d;
	      }
        float[] vs = find_minimum(ic,ie-ib-1,ib+1,c,what[i2][i1],prob);
	      next[i1]= vs[0];
        what[i2][i1] = vs[1];
	    }
	    for (int i1=0; i1<n1; i1++) {
	      prev[i1]=next[i1];
	      tt[i2][i1]=prev[i1];
	    }
    }
    backwardTrack(p, next, what);
    for (int i2=0; i2<n2;++i2) {
    for (int i1=0; i1<n1;++i1) {
      //tx[i2][i1]  = tt[i2][i1];
      tx[i1][i2]  = tt[i2][i1];
    }}
    return p;
  }



  private float[] find_minimum(
    int ic, int nc, int jc, float c, float pick, float[] prob)
  {
    float fm, f0, fp, a, b;
    if (0==ic) {
	    ic++;
	    fm=c;
	    f0=prob[ic];
	    fp=prob[ic+1];
    } else if (nc-1==ic) {
	    ic--;
	    fm=prob[ic-1];
	    f0=prob[ic];
	    fp=c;
    } else {
	    fm=prob[ic-1];
	    f0=c;
	    fp=prob[ic+1];
    }

    ic += jc;
    a = fm+fp-2f*f0;
    if (a <= 0.) { /* no minimum */
	    if (fm < f0 && fm < fp) {
	      pick = ic-1;
	      return new float[]{fm,pick};
	    } 
	    if (fp < f0 && fp < fm) {
	      pick = ic+1;
	      return new float[]{fp,pick};
	    } 
	    pick = ic;
	    return new float[]{f0,pick};
    }

    b = 0.5f*(fm-fp);
    a = b/a;
    if (a > 1.) {
	    pick = ic+1;
	    return new float[]{fp,pick};
    }

    if (a < -1.) {
	    pick = ic-1;
	    return new float[]{fm,pick};
    }

    if (f0 < 0.5*b*a) {
	    pick = ic;
	    return new float[]{f0,pick};
    }

    f0 -= 0.5*b*a;
    pick=ic+a;
    return new float[]{f0,pick};
  }

  private void forwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=0; i2<n2; i2++) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }


  private void backwardTrack(
    float[] path, float[] next, float[][] what)
  {
    float c, d, fc;
    int n1 = next.length;
    int n2 = path.length;
    c = FLT_MAX;
    fc = 0;
    /* minimum at the bottom */
    for (int i1=0; i1 < n1; i1++) {
	    d = next[i1];
	    if (d < c) {
	      c = d;
	      fc = i1;
	    }
    }
    /* coming up */
    for (int i2=n2-1; i2 >= 0; i2--) {
	    path[i2]=fc;
	    fc = interpolate(fc,i2,what);
    }
  }

  private float interpolate(float fc, int i2, float[][] what) {
    int n1 = what[0].length;
    int ic = round(fc-0.5f);
    fc -= ic;
    if (n1-1 <= ic) return what[i2][n1-1];
    if (0 > ic) return what[i2][0];
    fc = what[i2][ic]*(1f-fc)+what[i2][ic+1]*fc;
    return fc;
  }

  private void makeRhsWeightsInline(
    float[][] p23, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      float w23i = vel[i3][i2][k23];
      w23i *= w23i;
      ws[i3][i2] = w23i;
      b[i3][i2] = p23[i2][i3]*w23i;
    }}
  }


  private void makeRhsWeights(
    float[][] p23, float[][] p32, float[][][] vel, float[][] b, float[][] ws) 
  {
    int n3 = vel.length;
    int n2 = vel[0].length;
    int n1 = vel[0][0].length;
    for (int i3=0; i3<n3; ++i3) {
    for (int i2=0; i2<n2; ++i2) {
      int k23 = round(p23[i2][i3]);
      int k32 = round(p32[i3][i2]);
      k23 = max(0,k23);
      k23 = min(n1-1,k23);
      k32 = max(0,k32);
      k32 = min(n1-1,k32);
      float w23i = vel[i3][i2][k23];
      float w32i = vel[i3][i2][k32];
      w23i *= w23i;
      w32i *= w32i;
      ws[i3][i2] = w23i+w32i;
      b[i3][i2] = p23[i2][i3]*w23i+p32[i3][i2]*w32i;
    }}
  }


  private void makeRhsWeights(
    float[] p1, float[] p2, float[][] vel, float[] b, float[] ws) 
  {
    int n2 = vel.length;
    int n1 = vel[0].length;
    for (int i1=0; i1<n1; ++i1) {
      int k12 = round(p1[i1]);
      int k22 = round(p2[i1]);
      k12 = max(0,k12);
      k12 = min(n2-1,k12);
      k22 = max(0,k22);
      k22 = min(n2-1,k22);
      float w1i = vel[k12][i1];
      float w2i = vel[k22][i1];
      w1i *= w1i;
      w2i *= w2i;
      ws[i1] = w1i+w2i;
      b[i1] = p1[i1]*w1i+p2[i1]*w2i;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _gate;
  private float _an;
}
