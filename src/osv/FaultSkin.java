/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package osv;

import java.io.*;
import java.util.*;

import edu.mines.jtk.io.*;
import edu.mines.jtk.sgl.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A linked list of fault cells that may be used to analyze faults.
 *
 * @author Dave Hale and Xinming Wu, Colorado School of Mines
 * @version modified from Dave's FaultSkin in his ipf package
 */

public class FaultSkin implements Iterable<FaultCell>,Serializable {
  private static final long serialVersionUID = 1L;

  /**
   * Gets the cell that was the seed used to grow this skin.
   * @return the seed cell.
   */
  public FaultCell getSeed() {
    return _seed;
  }

  /**
   * Returns the number of cells in this skin.
   */
  public int size() {
    return _cellList.size();
  }

  /**
   * Gets an array of cells in this skin.
   * @return array of cells.
   */
  public FaultCell[] getCells() {
    return _cellList.toArray(new FaultCell[0]);
  }

  /**
   * Gets all cells in the specified skins.
   * @param skins array of skins for which to get cells.
   * @return array of cells.
   */
  public static FaultCell[] getCells(FaultSkin[] skins) {
    int ncell = countCells(skins);
    FaultCell[] cells = new FaultCell[ncell];
    int icell = 0;
    for (FaultSkin skin:skins)
      for (FaultCell cell:skin)
        cells[icell++] = cell;
    return cells;
  }

  /**
   * Gets a thinned fault likelihood image from fault skins.
   * @param skins array of skins for which to get a fault likelihood image.
   * @return array of a thinned fault likelihood image.
   */

  public static float[][][] getFl(int n1, int n2, int n3, FaultSkin[] skins) {
    float[][][] fl = new float[n3][n2][n1];
    for (FaultSkin skin:skins) {
    for (FaultCell cell:skin) {
      int i1 = cell.getI1();
      int i2 = cell.getI2();
      int i3 = cell.getI3();
      fl[i3][i2][i1] = cell.getFl();
    }}
    return fl;
  }

  /**
   * Returns the total number of cells in the specified skins.
   * @param skins array of skins for which to count cells.
   * @return the total number of cells.
   */
  public static int countCells(FaultSkin[] skins) {
    int ncell = 0;
    for (FaultSkin skin:skins)
      ncell += skin.size();
    return ncell;
  }

  /**
   * Returns an iterator for the cells in this skin. 
   * @return cell iterator.
   */
  public Iterator<FaultCell> iterator() {
    return _cellList.iterator();
  }

  /**
   * Returns array of arrays of cells linked above and below.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultCell[][] getCellsAB() {
    if (_cellsAB!=null)
      return _cellsAB;

    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search above for the top cell.
        FaultCell c = cell;
        for (FaultCell ca=c.ca; ca!=null; ca=c.ca)
          c = ca;

        // Add the top cell and all cells below it to a new list.
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        for (; c!=null; c=c.cb) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsAB = cellsList.toArray(new FaultCell[0][]);
    return _cellsAB;
  }

  /**
   * Returns array of arrays of cells linked left and right.
   * @return array of arrays of linked cells; by reference, not by copy.
   */
  public FaultCell[][] getCellsLR() {
    if (_cellsLR!=null)
      return _cellsLR;

    HashSet<FaultCell> cellSet = new HashSet<FaultCell>(size());
    ArrayList<FaultCell[]> cellsList = new ArrayList<FaultCell[]>();

    // For all cells in this skin, ...
    for (FaultCell cell:_cellList) {

      // If the cell is not already in an array, ...
      if (!cellSet.contains(cell)) {

        // Search left until we find no left nabor or until that left
        // nabor is the cell with which we began the search.
        FaultCell c = cell;
        for (FaultCell cl=c.cl; cl!=null && cl!=cell; cl=c.cl)
          c = cl;

        // Remember the leftmost cell found and add it to a new list.
        FaultCell cLeft = c;
        ArrayList<FaultCell> cList = new ArrayList<FaultCell>();
        cList.add(c);
        cellSet.add(c);

        // Add cells found to the right. Again beware of cycles.
        for (c=c.cr; c!=null && c!=cLeft; c=c.cr) {
          cList.add(c);
          cellSet.add(c);
        }

        // Convert the list to an array and add it to the list of arrays.
        cellsList.add(cList.toArray(new FaultCell[0]));
      }
    }
    assert _cellList.size()==cellSet.size();

    // Convert the list of arrays to the array of arrays to be returned.
    _cellsLR = cellsList.toArray(new FaultCell[0][]);
    checkCellArrays(_cellsLR);
    return _cellsLR;
  }

  /**
   * Smooths the normal vectors of cells in this skin.
   * @param nsmooth the number of smoothings.
   */
  public void smoothCellNormals(int nsmooth) {
    FaultCell.GetN getter = new FaultCell.GetN() {
      public float[] get(FaultCell cell) {
        return new float[]{cell.w1,cell.w2,cell.w3};
      }
    };
    FaultCell.SetN setter = new FaultCell.SetN() {
      public void set(FaultCell cell, float[] w) {
        float w1 = w[0]; 
        float w2 = w[1]; 
        float w3 = w[2];
        float ws = 1.0f/sqrt(w1*w1+w2*w2+w3*w3);
        w1 *= ws;
        w2 *= ws;
        w3 *= ws;
        cell.setNormalVector(w1,w2,w3);
      }
    };
    for (int ismooth=0; ismooth<nsmooth; ++ismooth)
      smoothN(getter,setter);
  }

  /**
   * Gets a cell nearest the centroid of this skin.
   * In illustrations, this cell is often a good representative.
   * @return the cell nearest the centroid.
   */
  public FaultCell getCellNearestCentroid() {
    float c1 = 0.0f;
    float c2 = 0.0f;
    float c3 = 0.0f;
    float cs = 0.0f;
    for (FaultCell c:_cellList) {
      c1 += c.fl*c.x1;
      c2 += c.fl*c.x2;
      c3 += c.fl*c.x3;
      cs += c.fl;
    }
    c1 /= cs;
    c2 /= cs;
    c3 /= cs;
    float dmin = Float.MAX_VALUE;
    FaultCell cmin = null;
    for (FaultCell c:_cellList) {
      float d = c.distanceSquaredTo(c1,c2,c3);
      if (d<dmin) {
        cmin = c;
        dmin = d;
      }
    }
    return cmin;
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault likelihoods.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   * @param lhc true, if left-handed coordinate system; false, otherwise.
   */
  public float[][] getCellXyzUvwRgbForLikelihood(
      float size, ColorMap cmap, boolean lhc) {
    return FaultCell.getXyzUvwRgbForLikelihood(size,cmap,getCells(),lhc);
  }

  /**
   * Gets arrays {xyz,uvw,rgb} of cell coordinates, normals and colors. In
   * these arrays, cells in this skin are represented by quads with specified
   * size, and colors corresponding to fault throws.
   * @param size size (in samples) of the quads.
   * @param cmap colormap used to compute rgb colors from cell properties.
   * @param lhc true, if left-handed coordinate system; false, otherwise.
   */
  public float[][] getCellXyzUvwRgbForThrow(
      float size, ColorMap cmap, boolean lhc) {
    return FaultCell.getXyzUvwRgbForThrow(size,cmap,getCells(),lhc);
  }

  public TriangleGroup getTriMesh(ColorMap cmap) {
    int nc = size();
    float[] rgb = new float[nc*6];
    float[] xyz = new float[nc*2*9];
    int k = 0;
    int i = 0;
    for (FaultCell cell:_cellList) {
      FaultCell cr = cell.cr;
      FaultCell cb = cell.cb;
      FaultCell rb = null;
      FaultCell br = null;
      if(cr!=null) rb = cr.cb;
      if(cb!=null) br = cb.cr;
      if(cr!=null&&rb!=null) {
        xyz[k++] = cell.x3;
        xyz[k++] = cell.x2;
        xyz[k++] = cell.x1;
        rgb[i++] = cell.fl;

        xyz[k++] = rb.x3;
        xyz[k++] = rb.x2;
        xyz[k++] = rb.x1;
        rgb[i++] = rb.fl;

        xyz[k++] = cr.x3;
        xyz[k++] = cr.x2;
        xyz[k++] = cr.x1;
        rgb[i++] = cr.fl;
      }
      if(cb!=null&&br!=null) {
        xyz[k++] = cell.x3;
        xyz[k++] = cell.x2;
        xyz[k++] = cell.x1;
        rgb[i++] = cell.fl;

        xyz[k++] = cb.x3;
        xyz[k++] = cb.x2;
        xyz[k++] = cb.x1;
        rgb[i++] = cb.fl;

        xyz[k++] = br.x3;
        xyz[k++] = br.x2;
        xyz[k++] = br.x1;
        rgb[i++] = br.fl;
      }
    }
    xyz = copy(k,0,xyz);
    rgb = copy(i,0,rgb);
    return new TriangleGroup(true,xyz,cmap.getRgbFloats(rgb));
  }

  public QuadGroup getQuadMeshStrike(ColorMap cmap) {
    int nc = size();
    float[] rgb = new float[nc*6];
    float[] xyz = new float[nc*2*9];
    int k = 0;
    int i = 0;
    for (FaultCell cell:_cellList) {
      FaultCell cb = cell.cb;
      FaultCell cr = cell.cr;
      FaultCell br = null;
      if(cb!=null) br = cb.cr;
      float fp = cell.fp;
      if(cr!=null&&br!=null&&cr!=null) {
        xyz[k++] = cell.x3;
        xyz[k++] = cell.x2;
        xyz[k++] = cell.x1;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = cb.x3;
        xyz[k++] = cb.x2;
        xyz[k++] = cb.x1;
        fp = cb.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = br.x3;
        xyz[k++] = br.x2;
        xyz[k++] = br.x1;
        fp = br.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = cr.x3;
        xyz[k++] = cr.x2;
        xyz[k++] = cr.x1;
        fp = cr.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;
      }
    }
    xyz = copy(k,0,xyz);
    rgb = copy(i,0,rgb);
    return new QuadGroup(true,xyz,cmap.getRgbFloats(rgb));
  }

  public TriangleGroup getTriMeshStrike(ColorMap cmap) {
    int nc = size();
    float[] rgb = new float[nc*6];
    float[] xyz = new float[nc*2*9];
    int k = 0;
    int i = 0;
    for (FaultCell cell:_cellList) {
      FaultCell cr = cell.cr;
      FaultCell cb = cell.cb;
      FaultCell rb = null;
      FaultCell br = null;
      if(cr!=null) rb = cr.cb;
      if(cb!=null) br = cb.cr;
      float fp = cell.fp;
      if(cr!=null&&rb!=null) {
        xyz[k++] = cell.x3;
        xyz[k++] = cell.x2;
        xyz[k++] = cell.x1;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = rb.x3;
        xyz[k++] = rb.x2;
        xyz[k++] = rb.x1;
        fp = rb.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = cr.x3;
        xyz[k++] = cr.x2;
        xyz[k++] = cr.x1;
        fp = cr.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;
      }
      if(cb!=null&&br!=null) {
        xyz[k++] = cell.x3;
        xyz[k++] = cell.x2;
        xyz[k++] = cell.x1;
        fp = cell.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = cb.x3;
        xyz[k++] = cb.x2;
        xyz[k++] = cb.x1;
        fp = cb.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;

        xyz[k++] = br.x3;
        xyz[k++] = br.x2;
        xyz[k++] = br.x1;
        fp = br.fp;
        if(fp>180) fp = 360-fp;
        rgb[i++] = fp;
      }
    }
    xyz = copy(k,0,xyz);
    rgb = copy(i,0,rgb);
    return new TriangleGroup(true,xyz,cmap.getRgbFloats(rgb));
  }

  public void updateStrike() {
    for (FaultCell cell:_cellList) {
      FaultCell cl = cell.cl;
      FaultCell cr = cell.cr;
      float fpl = 0.0f;
      float fpr = 0.0f;
      float scs = 0.0f;
      if(cl!=null) {
        float d2 = cell.x2-cl.x2;
        float d3 = cell.x3-cl.x3;
        float ds = sqrt(d2*d2+d3*d3);
        fpl = (float)Math.toDegrees(acos(d3/ds));
        scs += 1f;
      }
      if(cr!=null) {
        float d2 = -cell.x2+cr.x2;
        float d3 = -cell.x3+cr.x3;
        float ds = sqrt(d2*d2+d3*d3);
        fpr = (float)Math.toDegrees(acos(d3/ds));
        scs += 1f;
      }
      if(scs>0f) {cell.fp=(fpl+fpr)/scs;}
    }
  }

  public float getX1max() {
    float x1max = 0;
    for (FaultCell cell:_cellList) {
      if(cell.x1>x1max) x1max = cell.x1;
    }
    return x1max;
  }

  public float getX2max() {
    float x2max = 0;
    for (FaultCell cell:_cellList) {
      if(cell.x2>x2max) x2max = cell.x2;
    }
    return x2max;
  }

  public float getX3max() {
    float x3max = 0;
    for (FaultCell cell:_cellList) {
      if(cell.x3>x3max) x3max = cell.x3;
    }
    return x3max;
  }


  public void smooth(int n) { 
    for(int i=0; i<n; ++i) {
      for (FaultCell cell:_cellList) {
        FaultCell ca = cell.ca;
        FaultCell cb = cell.cb;
        FaultCell cl = cell.cl;
        FaultCell cr = cell.cr;
        float cs = 1f;
        float x1 = cell.x1;
        float x2 = cell.x2;
        float x3 = cell.x3;
        float fl = cell.fl;
        if(ca!=null) {
          cs += 1f;
          x1 += ca.x1;
          x2 += ca.x2;
          x3 += ca.x3;
          fl += ca.fl;
        }
        if(cb!=null) {
          cs += 1f;
          x1 += cb.x1;
          x2 += cb.x2;
          x3 += cb.x3;
          fl += cb.fl;
        }
        if(cl!=null) {
          cs += 1f;
          x1 += cl.x1;
          x2 += cl.x2;
          x3 += cl.x3;
          fl += cl.fl;
        }
        if(cr!=null) {
          cs += 1f;
          x1 += cr.x1;
          x2 += cr.x2;
          x3 += cr.x3;
          fl += cr.fl;
        }
        cs = 1f/cs;
        cell.x1 = x1*cs;
        cell.x2 = x2*cs;
        cell.x3 = x3*cs;
        cell.fl = fl*cs;
      }
    }
  }

  /**
   * Gets arrays of packed cell coordinates for cell links.
   * Each returned array contains packed (x,y,z) cell coordinates for
   * exactly one above-below or left-right linked list of cells.
   * @return array of arrays of packed xyz cell coordinates.
   */
  public float[][] getCellLinksXyz() {
    FaultCell[][] cellsAB = getCellsAB();
    FaultCell[][] cellsLR = getCellsLR();
    int nsAB = cellsAB.length; // number of segments for AB links
    int nsLR = cellsLR.length; // number of segments for LR links
    int ns = nsAB+nsLR; // total number of segments
    float[][] xyz = new float[ns][];
    float[][] rgb = new float[ns][];
    for (int is=0; is<ns; ++is) { // for all segments, ...
      FaultCell[] cells = (is<nsAB)?cellsAB[is]:cellsLR[is-nsAB]; // the cells
      int ncell = cells.length; // number of cells in this segment
      int np = ncell; // number of points in this segment
      if (is>=nsAB && cells[0].cl==cells[ncell-1]) // if a LR cycle, ...
        ++np; // then add one more so we end with the starting point
      float[] xyzi = new float[3*np]; // xyz for this segment
      for (int ip=0,ic=0; ip<np; ++ip) {
        FaultCell cell = cells[ip%ncell]; // % to handle any LR cycle
        xyzi[ic++] = cell.x3;
        xyzi[ic++] = cell.x2;
        xyzi[ic++] = cell.x1;
      }
      xyz[is] = xyzi;
    }
    return xyz;
  }

  /**
   * Returns a fault skin read from a file with specified name.
   * @param fileName the fault skin file name.
   * @return the fault skin.
   */
  public static FaultSkin readFromFile(String fileName) {
    FaultSkin skin = new FaultSkin();
    try {
      ArrayInputStream ais = new ArrayInputStream(fileName);
      int ncell = ais.readInt();
      ArrayList<FaultCell> cellList = new ArrayList<FaultCell>(ncell);
      for (int icell=0; icell<ncell; ++icell) {
        float x1 = ais.readFloat();
        float x2 = ais.readFloat();
        float x3 = ais.readFloat();
        float fl = ais.readFloat();
        float fp = ais.readFloat();
        float ft = ais.readFloat();
        FaultCell cell = new FaultCell(x1,x2,x3,fl,fp,ft);
        cell.skin = skin;
        cellList.add(cell);
        cell.s1 = ais.readFloat();
        cell.s2 = ais.readFloat();
        cell.s3 = ais.readFloat();
      }
      for (FaultCell cell:cellList) {
        int ida = ais.readInt();
        int idb = ais.readInt();
        int idl = ais.readInt();
        int idr = ais.readInt();
        if (ida!=INULL) 
          cell.ca = cellList.get(ida);
        if (idb!=INULL) 
          cell.cb = cellList.get(idb);
        if (idl!=INULL) 
          cell.cl = cellList.get(idl);
        if (idr!=INULL) 
          cell.cr = cellList.get(idr);
      }
      ais.close();
      skin._cellList = cellList;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return skin;
  }

  /**
   * Writes a fault skin to a file with specified name.
   * @param fileName the fault skin file name.
   * @param skin the fault skin.
   */
  public static void writeToFile(String fileName, FaultSkin skin) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fileName);
      aos.writeInt(skin.size());
      for (FaultCell cell:skin) {
        aos.writeFloat(cell.x1); 
        aos.writeFloat(cell.x2); 
        aos.writeFloat(cell.x3);
        aos.writeFloat(cell.fl); 
        aos.writeFloat(cell.fp); 
        aos.writeFloat(cell.ft);
        aos.writeFloat(cell.s1); 
        aos.writeFloat(cell.s2); 
        aos.writeFloat(cell.s3);
      }
      for (FaultCell cell:skin) {
        FaultCell[] nabors = new FaultCell[]{cell.ca,cell.cb,cell.cl,cell.cr};
        for (FaultCell nabor:nabors) {
          if (nabor!=null) {
            aos.writeInt(nabor.id);
          } else {
            aos.writeInt(INULL);
          }
        }
      }
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static FaultSkin readFromFileSlow(String fileName) {
    try {
      FileInputStream fis = new FileInputStream(fileName);
      ObjectInputStream ois = new ObjectInputStream(fis);
      FaultSkin skin = (FaultSkin)ois.readObject();
      ois.close();
      return skin;
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  public static void writeToFileSlow(String fileName, FaultSkin skin) {
    try {
      FileOutputStream fos = new FileOutputStream(fileName);
      ObjectOutputStream oos = new ObjectOutputStream(fos);
      oos.writeObject(skin);
      oos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // package

  /**
   * Constructs an empty skin.
   */
  FaultSkin() {
    _cellList = new ArrayList<FaultCell>();
  }

  /**
   * Adds the specified skinless cell to this skin.
   * @param cell the cell to be added.
   */
  void add(FaultCell cell) {
    assert cell.skin==null;
    cell.skin = this;
    if (_seed==null)
      _seed = cell;
    _cellList.add(cell);
    _cellsAB = null;
    _cellsLR = null;
  }

  /**
   * Smooths one value stored in the cells of this skin. The value smoothed is
   * that accessed by the specified getter and setter. Each smoothed value is
   * an average of the values in a cell and its cell nabors. 
   */
  void smooth1(FaultCell.Get1 getter, FaultCell.Set1 setter) {
    int ncell = size();
    float[] vals = new float[ncell];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float valCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float valNabor = getter.get(cellNabor);
          vals[icell] += valCell+valNabor;
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float vali = vals[icell]/(cnti>0.0f?cnti:1.0f);
      setter.set(cell,vali);
    }
  }

  /**
   * Smooths multiple values stored in the cells of this skin. The values
   * smoothed are those accessed by the specified getter and setter. Each
   * smoothed value is an average of the values in a cell and its cell nabors. 
   */
  void smoothN(FaultCell.GetN getter, FaultCell.SetN setter) {
    int ncell = size();
    int nval = getter.get(_seed).length;
    float[][] vals = new float[ncell][nval];
    float[] cnts = new float[ncell];
    FaultCell[] cellNabors = new FaultCell[4];
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float[] valsCell = getter.get(cell);
      cellNabors[0] = cell.ca;
      cellNabors[1] = cell.cb;
      cellNabors[2] = cell.cl;
      cellNabors[3] = cell.cr;
      for (FaultCell cellNabor:cellNabors) {
        if (cellNabor!=null) {
          float[] valsNabor = getter.get(cellNabor);
          for (int ival=0; ival<nval; ++ival)
            vals[icell][ival] += valsCell[ival]+valsNabor[ival];
          cnts[icell] += 2.0f;
        }
      }
    }
    for (int icell=0; icell<ncell; ++icell) {
      FaultCell cell = _cellList.get(icell);
      float cnti = cnts[icell];
      float scli = 1.0f/(cnti>0.0f?cnti:1.0f);
      for (int ival=0; ival<nval; ++ival)
        vals[icell][ival] *= scli;
      setter.set(cell,vals[icell]);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // private

  private static final int INULL = -Integer.MAX_VALUE; // null index

  private FaultCell _seed; // cell in this skin with highest fl; null, if empty
  private ArrayList<FaultCell> _cellList; // list of cells in this skin
  private FaultCell[][] _cellsAB; // arrays of cells from above to below
  private FaultCell[][] _cellsLR; // arrays of cells from left to right

  private void checkCellArrays() {
    if (_cellsAB!=null)
      checkCellArrays(_cellsAB);
    if (_cellsLR!=null)
      checkCellArrays(_cellsLR);
  }

  private static void checkCellArrays(FaultCell[][] cells) {
    HashSet<FaultCell> cellSet = new HashSet<FaultCell>();
    int ncell = cells.length;
    for (int icell=0; icell<ncell; ++icell) {
      int mcell = cells[icell].length;
      for (int jcell=0; jcell<mcell; ++jcell) {
        FaultCell c = cells[icell][jcell];
        assert cellSet.add(c);
      }
    }
  }

  private static void trace(String s) {
    System.out.println(s);
  }
}
