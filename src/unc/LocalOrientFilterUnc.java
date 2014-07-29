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

import util.*;
/**
 * compute two seismic normal vector fields from two structure-tensor fields
 * constructed with vertical causal and anti-causal filters
 * @author Xinming Wu and Dave Hale, Colorado School of Mines
 * @version 2014.04.12
 */
public class LocalOrientFilterUnc {
    
  /**
   * @param sigma1 half-width of window in vertical direction.
   * @param sigma2 half-width of window in direction smoothing along reflectors.
   * @param shift  shift step with which the input image is shift down or up
   * the image shifted down is used to construct structure tensors with vertical causal filter  
   * the image shifted up is used to construct structure tensors with vertical anticausal filter  
   */
  public LocalOrientFilterUnc(double sigma1, double sigma2, double shift) {
    _shift = (float)shift;
    _sigma1 = (float)sigma1;
    _lsfSmootherT = new LaplacianSmoother((int)sigma2);
    System.out.println("verticalSmoothing="+sigma1);
    System.out.println("literalSmoothing="+sigma2);
  }


  public void applyTest(float a, float[][] x, float[][] ys, float[][] ya) {
    causal(a,x,ys);
    anticausal(a,x,ya);
  }
  /**
   * For tests only
   */  
  public void applyForGradients(EigenTensors2 et, float[][] x,
    float[][] g1, float[][] g2, float[][] g11c, float[][] g12c,
    float[][] g22c, float[][] g11a, float[][] g12a, float[][] g22a) {
    int n2 = x.length;
    int n1 = x[0].length;
    Gradient.forward2D(x,g1,g2);
    // Gradient products.
    float[][] g11 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g22[i2][i1] = g2i*g2i;
        g12[i2][i1] = g1i*g2i;
      }
    }
    float[][][] gs = {g11,g22,g12};
    for (float[][] g:gs) {
      _lsfSmootherT.apply(_sd2,et,g,g);
    }
    causal(_sigma1,g11,g11c);
    causal(_sigma1,g12,g12c);
    causal(_sigma1,g22,g22c);
    anticausal2(_sigma1,g11c,g11c);
    anticausal2(_sigma1,g12c,g12c);
    anticausal2(_sigma1,g22c,g22c);
    anticausal(_sigma1,g11,g11a);
    anticausal(_sigma1,g12,g12a);
    anticausal(_sigma1,g22,g22a);
    causal2(_sigma1,g11a,g11a);
    causal2(_sigma1,g12a,g12a);
    causal2(_sigma1,g22a,g22a);

  }
  /**
   * Applies this filter to estimate normal vectors (1st eigenvectors).
   * @param x input array for 2-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   */
  public void applyForNormal(EigenTensors2 et, float[][] x,
    float[][] u1, float[][] u2, boolean causal)
  {
    // Gradient.
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] g1 = new float[n2][n1];
    float[][] g2 = new float[n2][n1];
    //float[][] xs = new float[n2][n1];
    //if(causal) southShift(_shift,x,xs); 
    //else       northShift(_shift,x,xs); 
    Gradient.forward2D(x,g1,g2);
    // Gradient products.
    float[][] g11 = new float[n2][n1];
    float[][] g12 = new float[n2][n1];
    float[][] g22 = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float g1i = g1[i2][i1];
        float g2i = g2[i2][i1];
        g11[i2][i1] = g1i*g1i;
        g22[i2][i1] = g2i*g2i;
        g12[i2][i1] = g1i*g2i;
      }
    }
    // Smoothed gradient products comprise the structure tensor.
    float[][] h = new float[n2][n1];
    float[][][] gs = {g11,g22,g12};
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(1.0);
    for (float[][] g:gs) {
      //ref.apply2(g,g);
      _lsfSmootherT.apply(_sd2,et,g,g);
    }
    eigenDecomposition(g11,g12,g22,u1,u2,causal);
  }
  // Compute eigenvectors, eigenvalues, and outputs that depend on them.
  private void eigenDecomposition(float[][] g11,float[][] g12,
    float[][] g22,float[][] u1,float[][] u2,boolean causal) 
  {
    int n2 = g11.length;
    int n1 = g11[0].length;
    if(causal) {
      causal(_sigma1,g11,g11);
      causal(_sigma1,g12,g12);
      causal(_sigma1,g22,g22);
      //anticausal2(_sigma1,g11,g11);
      //anticausal2(_sigma1,g12,g12);
      //anticausal2(_sigma1,g22,g22);

    } else {
      anticausal(_sigma1,g11,g11);
      anticausal(_sigma1,g12,g12);
      anticausal(_sigma1,g22,g22);
      //causal2(_sigma1,g11,g11);
      //causal2(_sigma1,g12,g12);
      //causal2(_sigma1,g22,g22);
    }
    float[][] a = new float[2][2];
    float[][] z = new float[2][2];
    float[] e = new float[2];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        a[0][0] = g11[i2][i1];
        a[0][1] = g12[i2][i1];
        a[1][0] = g12[i2][i1];
        a[1][1] = g22[i2][i1];
        Eigen.solveSymmetric22(a,z,e);
        float u1i = z[0][0];
        float u2i = z[0][1];
        if (u1i<0.0f) {
          u1i = -u1i;
          u2i = -u2i;
        }
        u1[i2][i1] = u1i;
        u2[i2][i1] = u2i;
      }
    }
  }

  /**
   * Applies this filter to estimate normal vectors (1st eigenvectors).
   * @param x input array for 3-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   * @param u3 3rd component of normal vector.
   */
  public void applyForNormal(EigenTensors3 et, float[][][] x, 
    float[][][] u1, float[][][] u2, float[][][] u3, boolean causal) 
  {
    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    float[][][] xs = new float[n3][n2][n1];
    System.out.println("shift="+_shift);
    if (causal) {southShift(_shift,x,xs);}
    else        {northShift(_shift,x,xs);}
    Gradient.forward3D(xs,g1,g2,g3);
    // Gradient products.
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    // Smoothed gradient products comprise the structure tensor.
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    for (float[][][] g:gs) {
      _lsfSmootherT.apply(_sd3,et,g,g);
    }
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,u1,u2,u3,causal);
  }
  /**
   * Applies this filter to estimate normal vectors (1st eigenvectors).
   * @param x input array for 3-D image.
   * @param u1 1st component of normal vector.
   * @param u2 2nd component of normal vector.
   * @param u3 3rd component of normal vector.
   */
  public void applyForNormal(EigenTensors3 et, float[][][] x, 
    float[][][] uc1, float[][][] uc2, float[][][] uc3, 
    float[][][] ua1, float[][][] ua2, float[][][] ua3) 
  {
    // Gradient.
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    float[][][] g1 = new float[n3][n2][n1];
    float[][][] g2 = new float[n3][n2][n1];
    float[][][] g3 = new float[n3][n2][n1];
    RecursiveGaussianFilter rgf =  new RecursiveGaussianFilter(1.0);
    rgf.apply100(x,g1);
    rgf.apply010(x,g2);
    rgf.apply001(x,g3);
    //Gradient.backward3D(x,g1,g2,g3);
    // Gradient products.
    float[][][] g11 = new float[n3][n2][n1];
    float[][][] g22 = new float[n3][n2][n1];
    float[][][] g33 = new float[n3][n2][n1];
    float[][][] g12 = new float[n3][n2][n1];
    float[][][] g13 = new float[n3][n2][n1];
    float[][][] g23 = new float[n3][n2][n1];
    computeGradientProducts(g1,g2,g3,g11,g12,g13,g22,g23,g33);
    // Smoothed gradient products comprise the structure tensor.
    float[][][] h = new float[n3][n2][n1];
    float[][][][] gs = {g11,g22,g33,g12,g13,g23};
    for (float[][][] g:gs) {
      _lsfSmootherT.apply(_sd3,et,g,g);
    }
    // Compute eigenvectors, eigenvalues, and outputs that depend on them.
    solveEigenproblems(g11,g12,g13,g22,g23,g33,uc1,uc2,uc3,true);
    solveEigenproblems(g11,g12,g13,g22,g23,g33,ua1,ua2,ua3,false);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  
  // Vertically causal and anti-causal filters are implemented 
  // by one-side recursive exponential filters.  
  private void causal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for (int i1=0; i1<n1; ++i1){ 
      float yi = y[0][i1] = x[0][i1];
      for (int i2=1; i2<n2; ++i2) 
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }

  private void causal(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[0] = x[0];
    for (int i1=1; i1<n1; ++i1) 
      y[i1] = yi = a*yi + b*x[i1];
  }

  private void causal(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        causal(a,x[i2],y[i2]);
      }
    });
  }

  private void causal(final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        causal(a,x[i3],y[i3]);
      }
    });

  }

  private void anticausal(float a, float[] x, float[] y) {
    int n1 = x.length;
    float b = 1.0f - a;
    float yi = y[n1-1] = x[n1-1];
    for(int i1=n1-2; i1>=0; --i1)
      y[i1] = yi = a*yi + b*x[i1];
  }
  private void anticausal2(float a, float[][] x, float[][] y) {
    int n2 = x.length;
    int n1 = x[0].length;
    float b = 1.0f - a;
    for(int i1=0; i1<n1; ++i1) {
      float yi = y[n2-1][i1] = x[n2-1][i1];
      for(int i2=n2-2; i2>=0; --i2)
        y[i2][i1] = yi = a*yi + b*x[i2][i1];
    }
  }

  private void anticausal(final float a, final float[][] x, final float[][] y) {
    final int n2 = x.length;
    Parallel.loop(n2, new Parallel.LoopInt() {
      public void compute(int i2) {
        anticausal(a,x[i2],y[i2]);
      }
    });
  }

  private void anticausal(final float a, final float[][][] x, final float[][][] y) {
    final int n3 = x.length;
    Parallel.loop(n3, new Parallel.LoopInt() {
      public void compute(int i3) {
        anticausal(a,x[i3],y[i3]); 
      }
    });
  }
 
  private float _shift  = 0.0f;
  private float _sigma1 = 0.8f;
  private static float[][] _weight;
  private LaplacianSmoother _lsfSmootherT;
  private LaplacianSmoother.Direction2 _sd2=LaplacianSmoother.Direction2.V; 
  private LaplacianSmoother.Direction3 _sd3=LaplacianSmoother.Direction3.VW; 

  private void computeGradientProducts(
    final float[][][] g1, final float[][][] g2, final float[][][] g3,
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33)
  {
    final int n1 = g1[0][0].length;
    final int n2 = g1[0].length;
    final int n3 = g1.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] g1i = g1[i3][i2];
          float[] g2i = g2[i3][i2];
          float[] g3i = g3[i3][i2];
          float[] g11i = g11[i3][i2];
          float[] g12i = g12[i3][i2];
          float[] g13i = g13[i3][i2];
          float[] g22i = g22[i3][i2];
          float[] g23i = g23[i3][i2];
          float[] g33i = g33[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float g1ii = g1i[i1];
            float g2ii = g2i[i1];
            float g3ii = g3i[i1];
            g11i[i1] = g1ii*g1ii;
            g22i[i1] = g2ii*g2ii;
            g33i[i1] = g3ii*g3ii;
            g12i[i1] = g1ii*g2ii;
            g13i[i1] = g1ii*g3ii;
            g23i[i1] = g2ii*g3ii;
          }
        }
      }
    });
  }

  private void solveEigenproblems(
    float[][][] g11, float[][][] g12, float[][][] g13,
    float[][][] g22, float[][][] g23, float[][][] g33,
    float[][][] u1, float[][][] u2, float[][][] u3, boolean causal) 
  {
    int n3 = g11.length;
    int n2 = g11[0].length;
    int n1 = g11[0][0].length;
    if(causal) {
      causal(_sigma1,g11,g11);
      causal(_sigma1,g12,g12);
      causal(_sigma1,g13,g13);
      causal(_sigma1,g22,g22);
      causal(_sigma1,g23,g23);
      causal(_sigma1,g33,g33);
      applyParallel(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    } else {
      anticausal(_sigma1,g11,g11);
      anticausal(_sigma1,g12,g12);
      anticausal(_sigma1,g13,g13);
      anticausal(_sigma1,g22,g22);
      anticausal(_sigma1,g23,g23);
      anticausal(_sigma1,g33,g33);
      applyParallel(g11,g12,g13,g22,g23,g33,u1,u2,u3);
    }

  }
  private void applyParallel( 
    final float[][][] g11, final float[][][] g12, final float[][][] g13,
    final float[][][] g22, final float[][][] g23, final float[][][] g33,
    final float[][][] u1, final float[][][] u2, final float[][][] u3) 
  {
    final int n3 = g11.length;
    final int n2 = g11[0].length;
    final int n1 = g11[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        double[][] a = new double[3][3];
        double[][] z = new double[3][3];
        double[] e = new double[3];
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            a[0][0] = g11[i3][i2][i1];
            a[0][1] = g12[i3][i2][i1];
            a[0][2] = g13[i3][i2][i1];
            a[1][0] = g12[i3][i2][i1];
            a[1][1] = g22[i3][i2][i1];
            a[1][2] = g23[i3][i2][i1];
            a[2][0] = g13[i3][i2][i1];
            a[2][1] = g23[i3][i2][i1];
            a[2][2] = g33[i3][i2][i1];
            Eigen.solveSymmetric33(a,z,e);
            float u1i = (float)z[0][0];
            float u2i = (float)z[0][1];
            float u3i = (float)z[0][2];
            if (u1i<0.0f) {
              u1i = -u1i;
              u2i = -u2i;
              u3i = -u3i;
            }
            u1[i3][i2][i1] = u1i;
            u2[i3][i2][i1] = u2i;
            u3[i3][i2][i1] = u3i;
          }
        }
      }
    });
  }
  private void southShift(float s, float[] x, float[] y) {
    int n1 = x.length; 
    SincInterpolator si = new SincInterpolator();
    //si.setUniform(n1,1,0,x);
    for (int i1=0; i1<n1; ++i1) {
      float x1i = (float)i1-s;
      y[i1] = si.interpolate(n1,1,0,x,x1i); 
    }
  } 
  private void northShift(float s, float[] x, float[] y) {
    int n1 = x.length; 
    SincInterpolator si = new SincInterpolator();
    //si.setUniform(n1,1,0,x);
    for (int i1=0; i1<n1; ++i1) {
      float x1i = (float)i1+s;
      y[i1] = si.interpolate(n1,1,0,x,x1i); 
    }
  } 
  private void southShift(
    final float s, final float[][] x, final float[][] y)
  {
    final int n2 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        southShift(s,x[i2],y[i2]);
      }
    });
  }
  private void northShift(
    final float s, final float[][] x, final float[][] y)
  {
    final int n2 = x.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        northShift(s,x[i2],y[i2]);
      }
    });
  }
  private void southShift(
    final float s, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        southShift(s,x[i3],y[i3]);
      }
    });
  }
  private void northShift(
    final float s, final float[][][] x, final float[][][] y)
  {
    final int n3 = x.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        northShift(s,x[i3],y[i3]);
      }
    });
  }
  
}
