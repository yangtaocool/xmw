package util;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;

 // Xinming Wu, Colorado School of Mines
 // 2013.11.24

public class NearestGridder {

  // 2D version
  public NearestGridder (float[] f, float[] x1, float[] x2) {
    _f = f;
    _x1 = x1;
    _x2 = x2;
    float[][] x = new float[2][];
    x[0] = _x1;
    x[1] = _x2;
    _kdt = new KdTree(x);
  }
  // 3D version
  public NearestGridder (float[] f, float[] x1, float[] x2, float[] x3) {
    _f = f;
    _x1 = x1;
    _x2 = x2;
    _x3 = x3;
    float[][] x = new float[3][];
    x[0] = _x1;
    x[1] = _x2;
    x[2] = _x3;
    _kdt = new KdTree(x);
  }

  // 2D version
  public float[][] grid(
    final Sampling s1, final Sampling s2) 
  { 
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final float f1 = (float)s1.getFirst();
    final float f2 = (float)s2.getFirst();
    final float d1 = (float)s1.getDelta();
    final float d2 = (float)s2.getDelta();
    final float[][] y = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1i = f1+d1*i1; 
	float x2i = f2+d2*i2;
	float[] xq = new float[] {x1i,x2i};
        int ind = _kdt.findNearest(xq);
        y[i2][i1] = _f[ind];
      }
    }});
    return y;
  }

  // 3D version
  public float[][][] grid(
    final Sampling s1, final Sampling s2, final Sampling s3) 
  { 
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    final float f1 = (float)s1.getFirst();
    final float f2 = (float)s2.getFirst();
    final float f3 = (float)s3.getFirst();
    final float d1 = (float)s1.getDelta();
    final float d2 = (float)s2.getDelta();
    final float d3 = (float)s3.getDelta();
    final float[][][] y = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+d1*i1; 
	  float x2i = f2+d2*i2;
	  float x3i = f3+d3*i3;
	  float[] xq = new float[] {x1i,x2i,x3i};
          int ind = _kdt.findNearest(xq);
          y[i3][i2][i1] = _f[ind];
	}
      }
    }});
    return y;
  }

  // 2D version
  public float[][] getDistance(
    final Sampling s1, final Sampling s2) 
  { 
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final float f1 = (float)s1.getFirst();
    final float f2 = (float)s2.getFirst();
    final float d1 = (float)s1.getDelta();
    final float d2 = (float)s2.getDelta();
    final float[][] y = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1i = f1+d1*i1; 
	float x2i = f2+d2*i2;
	float[] xq = new float[] {x1i,x2i};
        int ind = _kdt.findNearest(xq);
	float x1p = _x1[ind];
	float x2p = _x2[ind];
	float d1i = x1i-x1p;
	float d2i = x2i-x2p;
        y[i2][i1] = (float)Math.sqrt(d1i*d1i+d2i*d2i);
      }
    }});
    return y;
  }

  // 3D version
  public float[][][] getDistance(
    final Sampling s1, final Sampling s2, final Sampling s3) 
  { 
    final int n1 = s1.getCount();
    final int n2 = s2.getCount();
    final int n3 = s3.getCount();
    final float f1 = (float)s1.getFirst();
    final float f2 = (float)s2.getFirst();
    final float f3 = (float)s3.getFirst();
    final float d1 = (float)s1.getDelta();
    final float d2 = (float)s2.getDelta();
    final float d3 = (float)s3.getDelta();
    final float[][][] y = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+d1*i1; 
	  float x2i = f2+d2*i2;
	  float x3i = f3+d3*i3;
	  float[] xq = new float[] {x1i,x2i,x3i};
          int ind = _kdt.findNearest(xq);
	  float x1p = _x1[ind];
	  float x2p = _x2[ind];
	  float x3p = _x3[ind];
	  float d1i = x1i-x1p;
	  float d2i = x2i-x2p;
	  float d3i = x3i-x3p;
          y[i3][i2][i1] = (float)Math.sqrt(d1i*d1i+d2i*d2i+d3i*d3i);
	}
      }
    }});
    return y;
  }
///////////////////////////////////////////////////////////////////////
// private

  private KdTree _kdt;
  private float[] _f;
  private float[] _x1;
  private float[] _x2;
  private float[] _x3;
}
