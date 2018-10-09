/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package osv;

import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Fault cells in a 3D sampling grid. Each grid sample indexed by (i1,i2,i3)
 * contains either one fault cell or null. The grid facilitates searches for
 * cell nabors in skins and fast iterations along fault traces tangent to
 * fault strike and fault curves tangent to fault dip.
 * <p> 
 * Grid indices need not (and typically do not) begin at zero. Index bounds
 * for a fault cell grid are determined by the minima and maxima of indices of
 * cells used to construct the grid.
 *
 * @author Dave Hale and Xinming Wu, Colorado School of Mines
 * @version modified from Dave's FaultCellGrid in his ipf package
 */
public class FaultCellGrid {

  /**
   * Constructs a fault grid for specified cells. Grid index bounds are
   * determined by the minimum and maximum indices of the specified cells.
   * @param cells array of cells to be included in the grid.
   */
  public FaultCellGrid(FaultCell[] cells) {
    int i1min = Integer.MAX_VALUE;
    int i2min = Integer.MAX_VALUE;
    int i3min = Integer.MAX_VALUE;
    int i1max = -i1min;
    int i2max = -i2min;
    int i3max = -i3min;
    for (FaultCell cell:cells) {
      if (cell.i1<i1min) i1min = cell.i1;
      if (cell.i2<i2min) i2min = cell.i2;
      if (cell.i3<i3min) i3min = cell.i3;
      if (cell.i1>i1max) i1max = cell.i1;
      if (cell.i2>i2max) i2max = cell.i2;
      if (cell.i3>i3max) i3max = cell.i3;
    }
    _j1 = i1min;
    _j2 = i2min;
    _j3 = i3min;
    _n1 = 1+i1max-i1min;
    _n2 = 1+i2max-i2min;
    _n3 = 1+i3max-i3min;
    _cells = new FaultCell[_n3][_n2][_n1];
    for (FaultCell cell:cells)
      set(cell);
  }


  /**
   * Constructs a fault grid with specified dimensions
   * @param n1 1st dimension of the grid
   * @param n2 2nd dimension of the grid
   * @param n3 3rd dimension of the grid
   */
  public FaultCellGrid(int n1, int n2, int n3) {
    _j1 = 0;
    _j2 = 0;
    _j3 = 0;
    _n1 = n1;
    _n2 = n2;
    _n3 = n3;
    _cells = new FaultCell[n3][n2][n1];
  }

  /**
   * Gets the number of cells in the 1st dimension.
   * @return the number of cells.
   */
  public int getN1() {
    return _n1;
  }

  /**
   * Gets the number of cells in the 2nd dimension.
   * @return the number of cells.
   */
  public int getN2() {
    return _n2;
  }

  /**
   * Gets the number of cells in the 3rd dimension.
   * @return the number of cells.
   */
  public int getN3() {
    return _n3;
  }

  /**
   * Gets the lower bound on grid indices in the 1st dimension.
   * @return the lower bound.
   */
  public int getI1Min() {
    return _j1;
  }

  /**
   * Gets the lower bound on grid indices in the 2nd dimension.
   * @return the lower bound.
   */
  public int getI2Min() {
    return _j2;
  }

  /**
   * Gets the lower bound on grid indices in the 3rd dimension.
   * @return the lower bound.
   */
  public int getI3Min() {
    return _j3;
  }

  /**
   * Gets the upper bound on grid indices in the 1st dimension.
   * @return the upper bound.
   */
  public int getI1Max() {
    return _j1+_n1-1;
  }

  /**
   * Gets the upper bound on grid indices in the 2nd dimension.
   * @return the upper bound.
   */
  public int getI2Max() {
    return _j2+_n2-1;
  }

  /**
   * Gets the upper bound on grid indices in the 3rd dimension.
   * @return the upper bound.
   */
  public int getI3Max() {
    return _j3+_n3-1;
  }

  /**
   * Gets the fault cell with specified indices, if any.
   * @param i1 sample index in 1st dimension.
   * @param i2 sample index in 2nd dimension.
   * @param i3 sample index in 3rd dimension.
   * @return the fault cell; null, if none or if indices are out of bounds.
   */
  public FaultCell get(int i1, int i2, int i3) {
    i1 -= _j1; 
    i2 -= _j2; 
    i3 -= _j3;
    if (0<=i1 && i1<_n1 && 
        0<=i2 && i2<_n2 && 
        0<=i3 && i3<_n3) {
      return _cells[i3][i2][i1];
    } else {
      return null;
    }
  }

  /**
   * Sets the specified fault cell. Uses the cell's {x1,x2,x3} coordinates to
   * determine the indices of the cell in this grid.
   * @param cell the fault cell.
   */
  public void set(FaultCell cell) {
    int i1 = cell.i1-_j1;
    int i2 = cell.i2-_j2;
    int i3 = cell.i3-_j3;
    _cells[i3][i2][i1] = cell;
  }

  /**
   * Sets fault cells within a fault skin. 
   * Uses the cell's {x1,x2,x3} coordinates to
   * determine the indices of the cell in this grid.
   * @param cell the fault cell.
   */
  public void set(FaultSkin skin) {
    for (FaultCell cell:skin)
      set(cell);
  }

  /**
   * Sets the specified fault cells. 
   * Uses the cell's {x1,x2,x3} coordinates to
   * determine the indices of the cell in this grid.
   * @param cell the fault cell.
   */
  public void set(FaultCell[] cells) {
    for (FaultCell cell:cells)
      set(cell);
  }

  /**
   * Marks the fault cells near the specified fault cell. 
   * Finds the nearby celles with a search box centered 
   * at the specified fault cell.
   * @param cell the specified fault cell.
   * @param d1 1st radius of the search box
   * @param d2 2nd radius of the search box
   * @param d3 3rd radius of the search box
   */
  public void setCellsInBox(
    FaultCell cell, int d1, int d2, int d3) {
    int c1 = cell.i1;
    int c2 = cell.i2;
    int c3 = cell.i3;
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) ci.used = true;
    }}}
  }

  /**
   * Finds the fault cell nearest to the center of a search box.
   * @param c1 1st coordinate of the center point
   * @param c2 2nd coordinate of the center point
   * @param c3 3rd coordinate of the center point
   * @param d1 1st radius of the search box
   * @param d2 2nd radius of the search box
   * @param d3 3rd radius of the search box
   */
  public FaultCell findCellInBox(
    int c1, int c2, int c3, int d1, int d2, int d3) {
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    float dx = 1000f;
    FaultCell cell = null;
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      float di = (c1-i1)*(c1-i1)+(c2-i2)*(c2-i2)+(c3-i3)*(c3-i3);
      if(ci.skin==null&& di<dx) {
        dx = di;
        cell = ci;
      }
    }}}
    return cell;
  }

  /**
   * Checks in the search box to see if there are fault cells 
   * with fault strike similar to the specified fault cell.
   * @param x1 1st coordinate of the specified fault cell
   * @param x2 2nd coordinate of the specified fault cell
   * @param x3 3rd coordinate of the specified fault cell
   * @param fp strike of the specified fault cell
   * @param d1 1st radius of the search box
   * @param d2 2nd radius of the search box
   * @param d3 3rd radius of the search box
   */
  public boolean findCellsInBox(
    float x1, float x2, float x3, float fp, int d1, int d2, int d3) {
    int c1 = round(x1);
    int c2 = round(x2);
    int c3 = round(x3);
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    if(fp>180) fp = 360-fp;
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) {
        float fpi = ci.fp;
        if(fpi>180) fpi = 360-fpi;
        float dp = abs(fp-fpi);
        if(dp<40f) return true;
      }
    }}}
    return false;
  }

  /**
   * Checks in the search box to see if there are fault cells 
   * with fault strike similar to the specified fault cell.
   * @param cell the specified fault cell
   * @param d1 1st radius of the search box
   * @param d2 2nd radius of the search box
   * @param d3 3rd radius of the search box
   */
  public boolean findCellsInBox(
    FaultCell cell, int d1, int d2, int d3) {
    int c1 = cell.i1;
    int c2 = cell.i2;
    int c3 = cell.i3;
    float fp = cell.fp;
    int b1 = max(c1-d1,0);
    int b2 = max(c2-d2,0);
    int b3 = max(c3-d3,0);
    int e1 = min(c1+d1,_n1-1);
    int e2 = min(c2+d2,_n2-1);
    int e3 = min(c3+d3,_n3-1);
    for (int i3=b3; i3<=e3; i3++) {
    for (int i2=b2; i2<=e2; i2++) {
    for (int i1=b1; i1<=e1; i1++) {
      FaultCell ci = _cells[i3][i2][i1];
      if(ci!=null) {
        float fpi = ci.fp;
        float dp = abs(fp-fpi);
        dp = min(dp,360-dp);
        if(dp<40f) return true;
      }
    }}}
    return false;
  }


  /**
   * Finds a fault cell above the specified cell. Searches for a cell above
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell above,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell above.
   * @return the cell above; null, if none.
   */
  public FaultCell findCellAbove(FaultCell cell) {
    if (cell==null) return null;
    if (cell.ca!=null) return cell.ca;
    return findCellAboveBelow(true,cell);
  }

  /**
   * Finds a fault cell below the specified cell. Searches for a cell below
   * that lies nearest to the line containing the specified cell and its dip
   * vector. If the specified cell is already linked to a nabor cell below,
   * this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell below.
   * @return the cell below; null, if none.
   */
  public FaultCell findCellBelow(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cb!=null) return cell.cb;
    return findCellAboveBelow(false,cell);
  }

  /**
   * Finds a fault cell left of the specified cell. Searches for a cell left
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * left, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell left.
   * @return the cell left; null, if none.
   */
  public FaultCell findCellLeft(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cl!=null) return cell.cl;
    return findCellLeftRight(true,cell);
  }

  /**
   * Finds a fault cell right of the specified cell. Searches for a cell right
   * that lies nearest to the line containing the specified cell and its
   * strike vector. If the specified cell is already linked to a nabor cell
   * right, this method skips the search and simply returns that nabor cell. 
   * @param cell the cell for which to find a cell right.
   * @return the cell right; null, if none.
   */
  public FaultCell findCellRight(FaultCell cell) {
    if (cell==null) return null;
    if (cell.cr!=null) return cell.cr;
    return findCellLeftRight(false,cell);
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _j1,_j2,_j3; // min cell indices
  private int _n1,_n2,_n3; // numbers of cells
  private FaultCell[][][] _cells; // array of cells

  private void init(FaultCell[] cells) {
    if (cells!=null) {
      for (FaultCell cell:cells) {
        int i1 = cell.i1;
        int i2 = cell.i2;
        int i3 = cell.i3;
        _cells[i3-_j3][i2-_j2][i1-_j1] = cell;
      }
    }
  }

  private FaultCell findCellAboveBelow(boolean above, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float u1 = cell.u1;
    float u2 = cell.u2;
    float u3 = cell.u3;
    int k1 = 1;
    if (above) {
      k1 = -k1;
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        FaultCell c = get(i1+k1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float du = d1*u1+d2*u2+d3*u3;
          if (du>0.0f) {
            d1 -= du*u1;
            d2 -= du*u2;
            d3 -= du*u3;
            float d = d1*d1+d2*d2+d3*d3; // squared distance to dip line
            if (d<dmin) {
              cmin = c;
              dmin = d;
            }
          }
        }
      }
    }
    return cmin;
  }

  private FaultCell findCellAboveBelowX(boolean above, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float u1 = cell.u1;
    float u2 = cell.u2;
    float u3 = cell.u3;
    //int k1 = 1;
    int[] k1s = new int[]{1,2,3};
    if (above) {
      k1s = new int[]{-1,-2,-3};
      u1 = -u1;
      u2 = -u2;
      u3 = -u3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int p1=0; p1<3; p1 ++) {
      int k1 = k1s[p1];
    for (int k3=-1; k3<=1; ++k3) {
      for (int k2=-1; k2<=1; ++k2) {
        FaultCell c = get(i1+k1,i2+k2,i3+k3);
        if (c!=null) {
          float d1 = c.x1-x1;
          float d2 = c.x2-x2;
          float d3 = c.x3-x3;
          float du = d1*u1+d2*u2+d3*u3;
          if (du>0.0f) {
            d1 -= du*u1;
            d2 -= du*u2;
            d3 -= du*u3;
            float d = d1*d1+d2*d2+d3*d3; // squared distance to dip line
            if (d<dmin) {
              cmin = c;
              dmin = d;
            }
          }
        }
      }
    }}
    return cmin;
  }


  // The search for a cell left or right is not so straightforward as for a
  // cell above or below. The specified cell has eight adjacent samples. We
  // want a cell that is both nearby and located in the strike direction (if
  // right) or opposite direction (if left). We therefore first look for the
  // best cell among the N, E, S, and W adjacent samples, because they are
  // likely to be nearest, and we do not want to skip over them. If and only
  // if we do not find any candidate cells located in the specified direction,
  // we then look among the NE, SE, SW, and NW samples.
  private static final int[] K2LR = { 0, 1, 0,-1, 1, 1,-1,-1};
  private static final int[] K3LR = { 1, 0,-1, 0, 1,-1,-1, 1};
  private FaultCell findCellLeftRightX(boolean left, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v1 = cell.v1;
    float v2 = cell.v2;
    float v3 = cell.v3;
    if (left) {
      v1 = -v1;
      v2 = -v2;
      v3 = -v3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int ik=0; ik<8; ++ik) {
      if (ik==4 && cmin!=null)
        break;
      int k2 = K2LR[ik];
      int k3 = K3LR[ik];
      FaultCell c = get(i1,i2+k2,i3+k3);
      if (c!=null) {
        float d1 = c.x1-x1;
        float d2 = c.x2-x2;
        float d3 = c.x3-x3;
        float dv = d1*v1+d2*v2+d3*v3;
        if (dv>0.0f) {
          d1 -= dv*v1;
          d2 -= dv*v2;
          d3 -= dv*v3;
          float d = d1*d1+d2*d2+d3*d3; // squared distance to strike line
          if (d<dmin) {
            cmin = c;
            dmin = d;
          }
        }
      }
    }
    return cmin;
  }

  private FaultCell findCellLeftRight(boolean left, FaultCell cell) {
    int i1 = cell.i1;
    int i2 = cell.i2;
    int i3 = cell.i3;
    float x1 = cell.x1;
    float x2 = cell.x2;
    float x3 = cell.x3;
    float v1 = cell.v1;
    float v2 = cell.v2;
    float v3 = cell.v3;
    if (left) {
      v1 = -v1;
      v2 = -v2;
      v3 = -v3;
    }
    FaultCell cmin = null;
    float dmin = Float.MAX_VALUE;
    for (int k3=-5; k3<=5; k3++) {
    for (int k2=-5; k2<=5; k2++) {
      if(k2==0&&k3==0){continue;}
      FaultCell c = get(i1,i2+k2,i3+k3);
      if (c!=null) {
        float d1 = c.x1-x1;
        float d2 = c.x2-x2;
        float d3 = c.x3-x3;
        float dv = d1*v1+d2*v2+d3*v3;
        if (dv>0.0f) {
          d1 -= dv*v1;
          d2 -= dv*v2;
          d3 -= dv*v3;
          float d = d1*d1+d2*d2+d3*d3; // squared distance to strike line
          if (d<dmin) {
            cmin = c;
            dmin = d;
          }
        }
      }
    }}
    return cmin;
  }

}
