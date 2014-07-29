package util;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;
// backward and forward finite difference for gradient computation
public class Gradient {

  public static void forward2D(float[][] x, float[][] g1, float[][] g2) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1-1; ++i1) {
        float xi = x[i2][i1  ];
        float xp = x[i2][i1+1];
        g1[i2][i1] = xp-xi;
      }
      g1[i2][n1-1] = g1[i2][n1-2];
    }
    for (int i2=0; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2  ][i1];
        float xp = x[i2+1][i1];
        g2[i2][i1] = xp-xi;
      }
    }
    copy(g2[n2-2],g2[n2-1]);
  }

  public static void backward2D(float[][] x, float[][] g1, float[][] g2) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float xi = x[i2][i1  ];
        float xm = x[i2][i1-1];
        g1[i2][i1] = xi-xm;
      }
      g1[i2][0] = g1[i2][1];
    }
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float xi = x[i2  ][i1];
        float xm = x[i2-1][i1];
        g2[i2][i1] = xi-xm;
      }
    }
    copy(g2[1],g2[0]);
  }

  public static void center2D(float[][] x, float[][] g1, float[][] g2) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1-1; ++i1) {
        float x1m = x[i2][i1-1];
        float x1p = x[i2][i1+1];
        g1[i2][i1] = (x1p-x1m)*0.5f;
      }
      g1[i2][0] = x[i2][1]-x[i2][0];
      g1[i2][n1-1] = x[i2][n1-1]-x[i2][n1-2];
    }
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x2m = x[i2-1][i1];
        float x2p = x[i2+1][i1];
        g2[i2][i1] = (x2p-x2m)*0.5f;
      }
    }
    sub(x[1],x[0],g2[0]);
    sub(x[n2-1],x[n2-2],g2[n2-1]);
  }
  public static void forward3D(final float[][][] x, 
    final float[][][] g1, final float[][][] g2, final float[][][] g3) 
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {forward2D(x[i3],g1[i3],g2[i3]);}
    }); 
    for (int i3=0; i3<n3-1; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float xi = x[i3 ][i2][i1];
          float xp = x[i3+1][i2][i1];
          g3[i3][i2][i1] = xp-xi;
        }
      }
    }
    copy(g3[n3-2],g3[n3-1]);
  }

  public static void backward3D(final float[][][] x, 
    final float[][][] g1, final float[][][] g2, final float[][][] g3) 
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        backward2D(x[i3],g1[i3],g2[i3]);
    }}); 
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float xi = x[i3 ][i2][i1];
          float xm = x[i3-1][i2][i1];
          g3[i3][i2][i1] = xi-xm;
        }
      }
    }
    copy(g3[1],g3[0]);
  }
  
  public static void center3D(final float[][][] x,
    final float[][][] g1, final float[][][] g2, final float[][][] g3)
  {
    int n3 = x.length;
    int n2 = x[0].length;
    int n1 = x[0][0].length;
    Parallel.loop(n3, new Parallel.LoopInt() {
      public void compute(int i3) {
        center2D(x[i3],g1[i3],g2[i3]);
    }});
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x3m = x[i3-1][i2][i1];
          float x3p = x[i3+1][i2][i1];
          g3[i3][i2][i1] = (x3p-x3m)*0.5f;
        }
      }
    }
    sub(x[1],x[0],g3[0]);
    sub(x[n3-1],x[n3-2],g3[n3-1]);
  }
}
