/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package unc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Structure-oriented smoothing filter, taken from LocalSemblanceFilter for JTK
 * @author Xinming Wu, Colorado School of Mines
 * @version 2014.04.12
 */
public class LaplacianSmoother {

  public LaplacianSmoother (double sigma) {
    _scale = (float)sigma*((float)sigma+1.0f)/6.0f;
  }
  public enum Direction2 {
    U,V,UV
  }
  public enum Direction3 {
    U,V,W,UV,UW,VW,UVW
  }
  
  public void apply(float[] f, float[] g) {
    if (_scale==0.0f) {
        copy(f,g);
    } else {
      _lsf.apply(_scale,f,g);
    }
  }
  public void apply(
    Direction2 d, EigenTensors2 t, float[][] f, float[][] g) 
  {
    if (_scale==0.0f) {
      copy(f,g);
    } else {
      int n1 = f[0].length;
      int n2 = f.length;
      float[][] au = new float[n2][n1];
      float[][] av = new float[n2][n1];
      float[][] sf = new float[n2][n1];
      t.getEigenvalues(au,av);
      setEigenvalues(d,t);
      _lsf.apply(t,_scale,f,sf);
      _lsf.applySmoothS(sf,g);
      t.setEigenvalues(au,av);
    }
  }
  public void apply(
    Direction3 d, EigenTensors3 t, float[][][] f, float[][][] g) 
  {
    if (_scale==0.0f) {
      copy(f,g);
    } else {
      int n1 = f[0][0].length;
      int n2 = f[0].length;
      int n3 = f.length;
      float[][][] au = new float[n3][n2][n1];
      float[][][] av = new float[n3][n2][n1];
      float[][][] aw = new float[n3][n2][n1];
      float[][][] sf = new float[n3][n2][n1];
      t.getEigenvalues(au,av,aw);
      setEigenvalues(d,t);
      _lsf.applySmoothL(_kmax,f,sf);
      //_lsf.applySmoothS(f,sf);
      _lsf.apply(t,_scale,sf,g);
      //_lsf.applySmoothS(sf,g);
      t.setEigenvalues(au,av,aw);
    }
  }
  private float _scale;
  private static final double _small = 0.01;
  private static final int _niter = 1300;
  private static final LocalDiffusionKernel _ldk = 
  //new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D22);
    new LocalDiffusionKernel(LocalDiffusionKernel.Stencil.D71);
  private static final LocalSmoothingFilter _lsf = 
    new LocalSmoothingFilter(_small,_niter,_ldk);
  private static final double _kmax = 0.35;


  private static void setEigenvalues(Direction2 d, EigenTensors2 t) {
    float au = 0.0001f;
    float av = 0.0001f;
    if (d==Direction2.U || d==Direction2.UV)
      au = 1.0f;
    if (d==Direction2.V || d==Direction2.UV)
      av = 1.0f;
    t.setEigenvalues(au,av);
  }
  private static void setEigenvalues(Direction3 d, EigenTensors3 t) {
    float au = 0.000f;
    float av = 0.000f;
    float aw = 0.000f;
    if (d==Direction3.U || 
        d==Direction3.UV || 
        d==Direction3.UW ||
        d==Direction3.UVW)
      au = 1.0f;
    if (d==Direction3.V || 
        d==Direction3.UV || 
        d==Direction3.VW ||
        d==Direction3.UVW)
      av = 1.0f;
    if (d==Direction3.W || 
        d==Direction3.UW || 
        d==Direction3.VW ||
        d==Direction3.UVW)
      aw = 1.0f;
    t.setEigenvalues(au,av,aw);
  }

  private static float[] like(float[] f) {
    return new float[f.length];
  }
  private static float[][] like(float[][] f) {
    return new float[f.length][f[0].length];
  }
  private static float[][][] like(float[][][] f) {
    return new float[f.length][f[0].length][f[0][0].length];
  }

  private static Direction2 orthogonal(Direction2 d) {
    if (d==Direction2.U)
      return Direction2.V;
    else
      return Direction2.U;
  }
  private static Direction3 orthogonal(Direction3 d) {
    if (d==Direction3.U)
      return Direction3.VW;
    else if (d==Direction3.V)
      return Direction3.UW;
    else if (d==Direction3.W)
      return Direction3.UV;
    else if (d==Direction3.UV)
      return Direction3.W;
    else if (d==Direction3.UW)
      return Direction3.V;
    else
      return Direction3.U;
  }
}
