/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package osv;

import static edu.mines.jtk.util.ArrayMath.*;
/**
 * A fault cell is an oriented point located on a fault. 
 *
 * @author Xinming Wu, University of Texas at Austin
 * @version 2017.07.15
 */

public class FaultCell2 {

  public FaultCell2(int i1,int i2, float fl, float fp) {
    _i1 = i1;
    _i2 = i2;
    _fl = fl;
    _fp = fp;
  }

  public int getI1() {
    return _i1;
  }
  public int getI2() {
    return _i2;
  }

  public int[] getIndex() {
    return new int[]{_i1,_i2,_i3};
  }

  public float getFl() {
    return _fl;
  }

  public float getFp() {
    return _fp;
  }

  public float[] getFaultNormal() {
    return faultNormalVectorFromStrike(_fp);
  }

  public float[] getFaultStrikeVector() {
    return faultStrikeVectorFromStrike(_fp);
  }


  /**
   * Returns fault strike vector for specified strike and dip angles.
   * The dip angle theta is not used, but is provided for consistency.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {v1,v2,v3} of components for strike vector.
   */
  public static float[] faultStrikeVectorFromStrike(double phi) {
    double p = toRadians(phi);
    double cp = cos(p);
    double sp = sin(p);
    float v1 = (float)-cp;
    float v2 = (float) sp;
    return new float[]{v1,v2};
  }

  /**
   * Returns fault normal vector for specified strike and dip angles.
   * @param phi fault strike angle, in degrees.
   * @param theta fault dip angle, in degrees.
   * @return array {w1,w2,w3} of components for normal vector.
   */
  public static float[] faultNormalVectorFromStrike(double phi) {
    double p = toRadians(phi);
    double cp = cos(p);
    double sp = sin(p);
    float u1 = (float)(sp);
    float u2 = (float)(cp);
    return new float[]{u1,u2};
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  int _i1,_i2,_i3; // cell indices
  float _fl,_fp,_ft; // likelihood, strike (phi) and dip (theta)

}
