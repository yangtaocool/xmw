/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package unc;

import he.*;
import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
  *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2013.01.20
 */
public class Unconformity {

  public void setForLofu(double sigma1, double sigma2) {
    _shift = 0.0;
    _lofuSigma1  = sigma1;
    _lofuSigma2  = sigma2;
    _lofu = new LocalOrientFilterUnc(_lofuSigma1,_lofuSigma2,_shift);
  }

  public void setForLofu(double sigma1, double sigma2, double shift) {
    _shift = shift;
    _lofuSigma1  = sigma1;
    _lofuSigma2  = sigma2;
    _lofu = new LocalOrientFilterUnc(_lofuSigma1,_lofuSigma2,_shift);
  }
  // for test only
  public void applyForNormals(EigenTensors2 et, float[][] f, float[][] uc1, 
    float[][] uc2, float[][] ua1, float[][] ua2) 
  {
    _lofu.applyForNormal(et,f,uc1,uc2,true); 
    _lofu.applyForNormal(et,f,ua1,ua2,false); 
  }
  // for test only
  public void applyForGradients(EigenTensors2 et, float[][] f,
    float[][] g1, float[][] g2, float[][] g11c, float[][] g12c,
    float[][] g22c, float[][] g11a, float[][] g12a, float[][] g22a) {
    _lofu.applyForGradients(et,f,g1,g2,g11c,g12c,g22c,g11a,g12a,g22a); 
  }
  // compute unconformity likelihoods
  public void uncLikelihood(EigenTensors2 et,float[][] f, float[][] u) {
    int n2 = u.length; 
    int n1 = u[0].length; 
    float[][] uc1 = new float[n2][n1];
    float[][] uc2 = new float[n2][n1];
    float[][] ua1 = new float[n2][n1];
    float[][] ua2 = new float[n2][n1];
    _lofu.applyForNormal(et,f,uc1,uc2,true); 
    _lofu.applyForNormal(et,f,ua1,ua2,false); 
    int nb = (int)_shift+5;
    float[][] bd = zerofloat(nb,n2);
    add(bd,1.0f,bd);
    mul(uc1,ua1,uc1);
    mul(uc2,ua2,uc2);
    add(uc1,uc2,u);
    //copy(nb,n2,0,0,bd,0,0,u);
    //copy(nb,n2,0,0,bd,n1-nb,0,u);
  }
  // compute unconformity likelihoods
  public void uncLikelihood(EigenTensors3 et, float[][][] f, float[][][] u) {
    int n3 = u.length; 
    int n2 = u[0].length; 
    int n1 = u[0][0].length; 
    float[][][] uc1 = new float[n3][n2][n1];
    float[][][] uc2 = new float[n3][n2][n1];
    float[][][] uc3 = new float[n3][n2][n1];
    float[][][] ua1 = new float[n3][n2][n1];
    float[][][] ua2 = new float[n3][n2][n1];
    float[][][] ua3 = new float[n3][n2][n1];
    if (_shift>0.0f) {
      _lofu.applyForNormal(et,f,uc1,uc2,uc3,true); 
      _lofu.applyForNormal(et,f,ua1,ua2,ua3,false); 
      int nb = (int)_shift+5;
      float[][][] bd = zerofloat(nb,n2,n3);
      add(bd,1.0f,bd);
      mul(uc1,ua1,uc1);
      mul(uc2,ua2,uc2);
      mul(uc3,ua3,uc3);
      add(uc1,uc2,u);
      add(uc3,u  ,u);
      copy(nb,n2,n3,0,0,0,bd,0,0,0,u);
      copy(nb,n2,n3,0,0,0,bd,n1-nb,0,0,u);
    } else {
      _lofu.applyForNormal(et,f,uc1,uc2,uc3,ua1,ua2,ua3); 
      mul(uc1,ua1,uc1);
      mul(uc2,ua2,uc2);
      mul(uc3,ua3,uc3);
      add(uc1,uc2,u);
      add(uc3,u  ,u);
    }
  }
  // find all ridges of unconformity likelihoods 
  public void thin(int sigma, float[][] u, float[][] ut1, float[][] ut2) {
    thinParallel(sigma,u,ut1,ut2);
    //extend(u,ut1);
    //extend(u,ut2);
  }
  public void thin(int sigma, float[][][] u, float[][][] ut1, float[][][] ut2) {
    thinParallel(sigma,u,ut1,ut2);
    //extend(u,ut2);
    //extend(u,ut1);
  }
  public void thin(int sigma, float[] u, float[] ut1, float[] ut2) {
    int n1 = u.length;
    for (int i1=sigma; i1<n1-sigma; ++i1) {
      float ui = u[i1];
      float sum1 = 0.0f;
      float sum2 = 0.0f;
      for (int i=1; i<=sigma; ++i) {
        float dumi = u[i1-i] - ui; 
        float dupi = u[i1+i] - ui; 
        sum1 += dumi+dupi;
        sum2 += abs(dumi)+abs(dupi);
      }
      if (sum1==sum2 && ui<0.95f) {
        ut1[i1] = ui;
        ut2[i1-1] = u[i1-1];
        ut2[i1+1] = u[i1+1];
        ut2[i1  ] = ui;
      }

      /*
      if (sum1==sum2 && ui<0.55f) {
        ut1[i1] = ui;
        ut2[i1-2] = u[i1-2];
        ut2[i1-1] = u[i1-1];
        ut2[i1+1] = u[i1+1];
        ut2[i1] = ui;
      }
      */
    }
  }

  // connect ajacent points to form unconformity surfaces for 3D images
  // x thinned unconformity likelihood with only 1 sample vertically on the ridges
  // y unconformity likelihoods
  // output: unconformity surfaces saved in surf[ns][n3][n2]
  // ns: number of surfaces 
  // surf[is]: depth of samples on a surface
  public float[][][] uncSurfer(float th, float[][][] x, float[][][] y) {
    float[][][] u = copy(x);
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    int is = 0;
    int[][] mark = new int[n3][n2];
    float[][][] surf = new float[100][n3][n2];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float ui = u[i3][i2][i1];  
          int ip = 0;
          if(ui>0.0f && ui<th) {
            u[i3][i2][i1] = 0.0f;
            mark[i3][i2] = i1;
            double ym = y[i3][i2][i1-1];
            double yi = y[i3][i2][i1  ];
            double yp = y[i3][i2][i1+1];
            surf[is][i3][i2] = findPeak(i1,ym,yi,yp);
            ip ++;
            while (sum(mark)>0) {
              int[] ind = findMark(mark);
              int i3t = ind[0]; 
              int i2t = ind[1]; 
              int i1t = ind[2]; 
              mark[i3t][i2t] = 0;
              ip = floodFill(ip,i1t,i2t,i3t,1,1,1,u,y,mark,surf[is]);
            }
            if(ip>=10) { 
              is = is + 1;
              System.out.println("ip="+ip);
              System.out.println("is="+is);
            }
          }
        }
      }
    }
    surf = copy(n2,n3,is,0,0,0,surf);
    int ns3 = surf.length;
    int ns2 = surf[0].length;
    int ns1 = surf[0][0].length;
    System.out.println("ns3="+ns3);
    System.out.println("ns2="+ns2);
    System.out.println("ns1="+ns1);
    return trianglesForSurface(y,surf);
  }

  private void thinParallel(
     final int sigma, final float[][] u, 
     final float[][] ut1, final float[][] ut2) 
  {
    final int n2 = u.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        thin(sigma,u[i2],ut1[i2],ut2[i2]);
      }
    });
  }

  private void thinParallel(
    final int sigma, final float[][][] u, 
    final float[][][] ut1, final float[][][] ut2) 
  {
    final int n3 = u.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        thin(sigma,u[i3],ut1[i3],ut2[i3]);
      }
    });
  }
  
  private void extend(float[][] u, float[][] ut) {
    int n2 = u.length;
    int n1 = u[0].length;
    int[][] mark = new int[n2][n1];
    for (int i2=1; i2<n2-1; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float um  = u[i2-1][i1];
        float up  = u[i2+1][i1];
        float utm = ut[i2+1][i1];
        float utp = ut[i2+1][i1];
        int marki = mark[i2][i1];
        float uti= ut[i2  ][i1];
        if(uti<1.0f && marki==0){
          if(utm<1.0f) {ut[i2-1][i1] = um;mark[i2-1][i1] = 1;}
          if(utp<1.0f) {ut[i2+1][i1] = up;mark[i2+1][i1] = 1;}
         }
      }
    }
  }
  
  private void extend(float[][][] u, float[][][] ut) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    for (int i3=0; i3<n3; ++i3) 
      extend(u[i3],ut[i3]);
    float[][] ui2  = new float[n3][n1];
    float[][] uti2 = new float[n3][n1];
    for (int i2=0; i2<n2; ++i2) { 
      for (int i3=0; i3<n3; ++i3) { 
        ui2[i3]  = u[i3][i2];
        uti2[i3] = ut[i3][i2];
      }
      extend(ui2,uti2);
      for (int i3=0; i3<n3; ++i3) 
        ut[i3][i2] = uti2[i3];
    }
  }
  
  private float[][][] trianglesForSurface(float[][][] u, float[][][] s) {
    int ns = s.length;
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    //SincInterpolator usi = new SincInterpolator();
    //usi.setUniform(n1,1,0,n2,1,0,n3,1,0,u);
    float[][][] sf = new float[2][ns][];
    float[] xyz = new float[n2*n3*18];
    float[] rgb = new float[n2*n3*6];
    for (int is=0; is<ns; ++is) {
      int[] k = new int[2];
      float[][] si = s[is];
      for (int i3=0; i3<n3-1; ++i3) {
        for (int i2=0; i2<n2-1; ++i2) {
          int[] sia = new int[4];
          int[] i3a = new int[4];
          int[] i2a = new int[4];
          int[] ind = new int[4];
          i3a[0] = i3;   i3a[1] = i3+1;
          i3a[2] = i3+1; i3a[3] = i3;
          i2a[0] = i2;   i2a[1] = i2;
          i2a[2] = i2+1; i2a[3] = i2+1;
          int nz = 0;
          for (int i=0; i<4; i++) {
            int i3t = i3a[i]; 
            int i2t = i2a[i]; 
            float sit = si[i3t][i2t];
            if (sit>0.f) {ind[nz] = i;nz++;}
          }
          if (nz<=2) continue;
          if (nz==3) 
            k = triangle(k,i2a,i3a,ind,u,si,xyz,rgb);
          if (nz==4) {
            k = triangle(k,i2a,i3a,ind,u,si,xyz,rgb);
            ind[1] = 2; ind[2] = 3;
            k = triangle(k,i2a,i3a,ind,u,si,xyz,rgb);
          }
        }
      }
      sf[0][is] = copy(k[0],0,xyz);
      ColorMap mp = new ColorMap(0.0,1.0,ColorMap.JET);
      rgb = sub(1.0f,rgb);
      sf[1][is] = mp.getRgbFloats(copy(k[1],0,rgb)); 
    } 
    return sf;
  } 

  private int[] triangle(
    int[] k, int[] i2a, int[] i3a, int[] ind, 
    float[][][] u, float[][] s, float[] xyz, float[] rgb)
  {
     int k1 = k[0];
     int k2 = k[1];
     int n3 = u.length;
     int n2 = u[0].length;
     int n1 = u[0][0].length;
     SincInterpolator usi = new SincInterpolator();
     for (int i=0; i<3; i++) {
       int idt = ind[i];
       int i3t = i3a[idt];
       int i2t = i2a[idt];
       float x1i = s[i3t][i2t];
       float x3i = (float) i3t;
       float x2i = (float) i2t;
       xyz[k1  ] =  x3i;
       xyz[k1+1] =  x2i;
       xyz[k1+2] =  x1i;
       rgb[k2] = usi.interpolate(n1,1,0,n2,1,0,n3,1,0,u,x1i,x2i,x3i);
       k1 = k1+3;
       k2 = k2+1;
     } 
     k[0] = k1;
     k[1] = k2;
     return k;
  }
  private int floodFill(
    int ip, int i1t, int i2t, int i3t, 
    int d1, int d2, int d3, 
    float[][][] u, float[][][] y, int[][] mark, float[][] surf) 
  {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    int i1b = i1t-d1; if(i1b<0) i1b=0;
    int i2b = i2t-d2; if(i2b<0) i2b=0;
    int i3b = i3t-d3; if(i3b<0) i3b=0;
    int i1e = i1t+d1; if(i1e>n1-1) i1e=n1-1;
    int i2e = i2t+d2; if(i2e>n2-1) i2e=n2-1;
    int i3e = i3t+d3; if(i3e>n3-1) i3e=n3-1;
    float d = 4.f;
    for (int i3=i3b; i3<=i3e; ++i3) { 
      for (int i2=i2b; i2<=i2e; ++i2) { 
        for (int i1=i1b; i1<=i1e; ++i1) { 
          float ui = u[i3][i2][i1];
          if(ui>0.0f && ui<0.85f) {
            float d1i = i1-i1t;
            float d2i = i2-i2t;
            float d3i = i3-i3t;
            float di  = d1i*d1i+d2i*d2i+d3i*d3i;
            if(di<d) {
              double ym = y[i3][i2][i1-1];
              double yi = y[i3][i2][i1  ];
              double yp = y[i3][i2][i1+1];
              surf[i3][i2] = findPeak(i1,ym,yi,yp);
              mark[i3][i2] = i1;
              u[i3][i2][i1] = 0.0f;
              ip++;
            }
          }
        }
      }
    }
    return ip;
  }
  private int[] findMark(int[][]mark) {
    int n3 = mark.length;
    int n2 = mark[0].length;
    int[] ind = new int[3];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        int i1 = mark[i3][i2];
        if (i1>0) {
          ind[0] = i3;
          ind[1] = i2;
          ind[2] = i1;
          break;
        }
      }
    }
    return ind;
  }

  private float findPeak(int i1, double u1, double u2, double u3) {
    float z = (float) i1;
    double a = u1-u3;
    double b = 2.0*(u1+u3)-4.0*u2;
    double d = a/b;
    return (z+(float)d);
  }
    
  private double _shift=0.0; 
  private double _lofuSigma1 = 0.95;
  private double _lofuSigma2 = 128.0;
  private LocalOrientFilterUnc _lofu = 
    new LocalOrientFilterUnc(_lofuSigma1,_lofuSigma2,_shift);
}
