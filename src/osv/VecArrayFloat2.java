/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package osv;

/**
 * A vector represented by a 2D array[n2][n1] of floats.
 * @author Dave Hale, Colorado School of Mines
 * @version 2009.09.15
 */
public class VecArrayFloat2 implements Vec {

  /**
   * Constructs a zero vector with specified dimensions.
   * @param n1 the number of floats in the 1st dimension.
   * @param n2 the number of floats in the 2nd dimension.
   */
  public VecArrayFloat2(int n1, int n2) {
    _a = new float[n2][n1];
    _n1 = n1;
    _n2 = n2;
  }

  /**
   * Constructs a vector that wraps the specified array of floats.
   * @param a the array of floats; by reference, not by copy.
   */
  public VecArrayFloat2(float[][] a) {
    _a = a;
    _n1 = a[0].length;
    _n2 = a.length;
  }

  /**
   * Gets the array of floats wrapped by this vector.
   * @return the array of floats; by reference, not by copy.
   */
  public float[][] getArray() {
    return _a;
  }

  /**
   * Gets the number of floats in the 1st array dimension.
   * @return the number of floats in the 1st dimension.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of floats in the 2nd array dimension.
   * @return the number of floats in the 2nd dimension.
   */
  public int getN2() {
    return _n2;
  }

  public double epsilon() {
    return Math.ulp(1.0f);
  }

  public VecArrayFloat2 clone() {
    VecArrayFloat2 v = new VecArrayFloat2(_n1,_n2);
    for (int i2=0; i2<_n2; ++i2)
      System.arraycopy(_a[i2],0,v._a[i2],0,_n1);
    return v;
  }

  public double dot(Vec vthat) {
    float[][] athis = _a;
    float[][] athat = ((VecArrayFloat2)vthat)._a;
    double sum = 0.0;
    for (int i2=0; i2<_n2; ++i2) {
      float[] athis2 = athis[i2];
      float[] athat2 = athat[i2];
      for (int i1=0; i1<_n1; ++i1)
        sum += athis2[i1]*athat2[i1];
    }
    return sum;
  }

  public double norm2() {
    double sum = 0.0;
    for (int i2=0; i2<_n2; ++i2) {
      float[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        double ai = a2[i1];
        sum += ai*ai;
      }
    }
    return Math.sqrt(sum);
  }

  public void zero() {
    for (int i2=0; i2<_n2; ++i2) {
      float[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        a2[i1] = 0.0f;
      }
    }
  }

  public void scale(double s) {
    for (int i2=0; i2<_n2; ++i2) {
      float[] a2 = _a[i2];
      for (int i1=0; i1<_n1; ++i1) {
        a2[i1] *= s;
      }
    }
  }

  public void add(double sthis, Vec vthat, double sthat) {
    float[][] athis = _a;
    float[][] athat = ((VecArrayFloat2)vthat)._a;
    float fthis = (float)sthis;
    float fthat = (float)sthat;
    for (int i2=0; i2<_n2; ++i2) {
      float[] athis2 = athis[i2];
      float[] athat2 = athat[i2];
      for (int i1=0; i1<_n1; ++i1) {
        athis2[i1] = athis2[i1]*fthis+athat2[i1]*fthat;
      }
    }
  }

  private float[][] _a;
  private int _n1,_n2;
}
