package hv;

import java.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Extract a single seismic horizon surface with control points
 * the constraints derived from control points are incorporated in a 
 * preconditioner in the conjugate gradient method used to solve the 
 * linear system for horizon extracting.
 * @author Xinming Wu and Dave Hale
 * @version 2014.03.12
 */

public class SurfaceExtractorC {
 
  // scale the curvature term 
  public void setWeights(float w){
    _weight = w;
  }
  
  public void setSmoothings(float sigma1, float sigma2){
    _sigma1 = sigma1;
    _sigma2 = sigma2;
  }
  
  public void setCG(float small, int niter){
    _small = small;
    _niter = niter;
  }
  
  public void setExternalIterations(int exniter){
    _exniter = exniter;
  }
  // find the peak or trough nearest to each control point
  public float[] refineConstraints(float[] k1, float[] k2, float[] k3, float[][][] u) {
    int np = k1.length;
    int n1 = u[0][0].length;
    for (int ip=0; ip<np; ++ip) {
      float k1i = k1[ip];
      int i1 = (int)k1i;
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      int i1m = i1-1; 
      int i1p = i1+1; 
      if(i1m<0  ){i1m=0;   i1=1;   i1p=2;}
      if(i1p>=n1){i1m=n1-3;i1=n1-2;i1p=n1-1;}
      float um = u[i3][i2][i1m];
      float ui = u[i3][i2][i1 ];
      float up = u[i3][i2][i1p];
      float kp = parabolicPeak(k1i,um,ui,up);
      if(abs(kp-k1i)<2.0 && k1i*kp>=0.0f){k1[ip]=kp;}
    }
    return k1;
  }  

  // Interpolate an initial surface passing through control points
  public float[][] surfaceInitialization
    (int n2, int n3, float lmt, float[] k1, float[] k2, float[] k3) 
  {
    if (k1.length==1) {
      float[][] surf = zerofloat(n2,n3);
      add(surf,k1[0],surf);
      return surf; 
    } else {
      Sampling s2 = new Sampling(n2,1.0f,0.0f);
      Sampling s3 = new Sampling(n3,1.0f,0.0f);
      RadialInterpolator2.Biharmonic bs = new RadialInterpolator2.Biharmonic();
      RadialGridder2 rg = new RadialGridder2(bs,k1,k2,k3);    
      float[][] surf = rg.grid(s2,s3);
      surfaceCorrect(surf,lmt);
      checkControlPoints(k2, k3, surf); 
      return surf;
    }
  }

  // Updates the surface using the seismic normal vectors and control points.
  public void surfaceUpdateFromSlopes
    (float[][][] ep, float[][][] p ,float[][][] q,
     float[] k1, float[] k2, float[] k3,float[][] surf)
  {	
    int n3 = p.length; 
    int n2 = p[0].length; 
    int n1 = p[0][0].length; 
    float lmt = (float)n1-1.f;
    float[][] surft = copy(surf);
    float avg = sum(surft)/(n2*n3);
    float[][] b   = new float[n3][n2]; 
    float[][] pi1 = new float[n3][n2]; 
    float[][] qi1 = new float[n3][n2]; 
    float[][] wi1 = new float[n3][n2]; 
    int np = k1.length;
    float[] cf = new float[np]; // force of constraints
    for (int n=1; n<=_exniter; n++){
      System.out.println(" Iteration "+n+"......");
      VecArrayFloat2 vb    = new VecArrayFloat2(b);
      VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
      updateSlopesAndWeights(p,q,ep,surf,pi1,qi1,wi1);
      System.out.println(" maxp="+max(pi1));
      System.out.println(" maxq="+max(qi1));
      System.out.println(" maxw="+max(wi1));
      A2 a2 = new A2(wi1,_weight);
      M2 m2 = new M2(avg,_sigma1,_sigma2,wi1,k2,k3);
      CgSolver cs = new CgSolver(_small, _niter);
      vb.zero();
      makeRhs(wi1,pi1,qi1,b);
      cs.solve(a2,m2,vb,vsurf);
      surf = vsurf.getArray();
      surfaceCorrect(surf,lmt);
      checkControlPoints(k2, k3, surf); 
      constraintForce(k2,k3,surft,surf,cf);
      float ad = sum(abs(sub(surft,surf)))/(n3*n2); 
      System.out.println(" Average adjustments per sample = "+ad);
      if (ad<0.02f) break;
      surft = copy(surf);
    }
  
    // show the constraint force for each control point,
    // the points with little constraint force are not important, 
    // and hence can be removed from the constraints.
    checkConstraintForce(k2,k3,cf);
  }

  public void surfaceRefine(float[][] surf, float[][][] u) {
    int n3 = u.length;
    int n2 = u[0].length;
    float[] k1 = new float[1];
    float[] k2 = new float[1];
    float[] k3 = new float[1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        k2[0] = (float)i2;
        k3[0] = (float)i3;
        k1[0] = surf[i3][i2];
        surf[i3][i2]=refineConstraints(k1,k2,k3,u)[0];
      }
    }
  }
  public float[][] updateWeights(float[][] surf, float[][][] w) {
    int n3 = w.length;
    int n2 = w[0].length;
    int n1 = w[0][0].length;
    float[][] wi = new float[n3][n2];
    SincInterpolator wsi = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
	      wi[i3][i2] = wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,w,x1i,x2i,x3i);
      }
    }
    return wi;
  }
  private static void updateSlopesAndWeights (
    float[][][] p, float[][][] q, float[][][] ep,
    float[][] surf, float[][] pi1, float[][] qi1,float[][] wi1)
  {
    int n3 = p.length;
    int n2 = p[0].length;
    int n1 = p[0][0].length;
    SincInterpolator psi = new SincInterpolator();
    SincInterpolator qsi = new SincInterpolator();
    SincInterpolator wsi = new SincInterpolator();
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        double x2i = (double)i2;
        double x3i = (double)i3;
        double x1i = (double)surf[i3][i2];
        float wi = wsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0,ep,x1i,x2i,x3i);
        wi1[i3][i2] = (wi>0.0005f)?wi:0.0005f;
        pi1[i3][i2] = psi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, p,x1i,x2i,x3i);
	      qi1[i3][i2] = qsi.interpolate(n1,1.0,0.0,n2,1.0,0.0,n3,1.0,0.0, q,x1i,x2i,x3i);
      }
    }
  }

  private void surfaceCorrect(float[][] surf, float lmt) {
    int n1 = surf[0].length;
    int n2 = surf.length;
    for(int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        if (surf[i2][i1]<0.f) surf[i2][i1]=0.f;
        if (surf[i2][i1]>lmt) surf[i2][i1]=lmt;
      }
    }
  }

  private static class A2 implements CgSolver.A{
    A2(float[][] wp, float w){
      _w  = w;
      _wp = wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      VecArrayFloat2 v2z = v2x.clone();
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      float[][] yt = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      applyLhs(_wp,x,y);
      if (_w>0.0f) {
        applyLhs(_wp,x,yt);
        applyLhs(_wp,yt,yy);
        v2y.add(1.f,v2yy,_w);
      }
    }
    private float[][] _wp;
    private float _w;
  }

   // Preconditioner; includes smoothers and (optional) constraints.
  private static class M2 implements CgSolver.A {
    M2(float avg, float sigma1, float sigma2, float[][] wp, float[] k2, float[] k3) {
      _avg = avg;
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _wp = wp;
      if (k2!=null && k3!=null) {
        _k2 = copy(k2);
        _k3 = copy(k3);
      }
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      copy(x,y);
      constrain(_k2,_k3,y);
      //removeAverage(_avg,y);
      smooth2(_sigma2,_wp,y);
      smooth1(2.f*_sigma1,_wp,y);
      smooth2(_sigma2,_wp,y);
      //removeAverage(_avg,y);
      constrain(_k2,_k3,y);
    }
    private float _avg,_sigma1,_sigma2;
    private float[][] _wp;
    private float[] _k2,_k3;
  }

  private static void constraintForce
    (float[] k2, float[] k3, float[][] fp, float[][] fi, float[] cf)
  {
    int np = k2.length; 
    int n3 = fp.length;
    int n2 = fp[0].length;
    for (int ip=0; ip<np; ++ip) {
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      int i2m = i2-1;
      int i2p = i2+1;
      int i3m = i3-1;
      int i3p = i3+1;
      float c  = 0.0f;
      float dw = 0.0f; 
      float de = 0.0f; 
      float ds = 0.0f; 
      float dn = 0.0f; 
      if(i3m>=0){c +=1.0f;dw = abs(fp[i3m][i2]-fi[i3m][i2]);}
      if(i3p<n3){c +=1.0f;de = abs(fp[i3p][i2]-fi[i3p][i2]);}
      if(i2m>=0){c +=1.0f;dn = abs(fp[i3][i2m]-fi[i3][i2m]);}
      if(i2p<n2){c +=1.0f;ds = abs(fp[i3][i2p]-fi[i3][i2p]);}
      cf[ip] = (dw+de+dn+ds)/c;
    }
  }
  private static void checkControlPoints(float[] k2, float[] k3, float[][] f) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip];
        int i3 = (int)k3[ip];
        System.out.println(" i2="+i2+" i3="+i3+" f1="+f[i3][i2]);
      }
    }
  }

  private static void removeAverage(float avg, float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    double n1n2 = (double)n1*(double)n2;
    double sumx = sum(x);
    double avgx = sumx/n1n2;
    float  avgi = (float)avgx-avg;
    for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) 
        x[i2][i1] -= avgi;
  }
  private static void constrain(float[] k2, float[] k3, float[][] x) {
    if (k2!=null && k3!=null) {
      int np = k2.length;
      for (int ip=0; ip<np; ++ip) {
        int i2 = (int)k2[ip]; 
        int i3 = (int)k3[ip]; 
        x[i3][i2] = 0.f;
      }
    }
  }

  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    float[] st = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        xt[i1] = x[i2][i1];
        st[i1] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] s, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        xt[i2] = x[i2][i1];
        st[i2] = s[i2][i1];
      }
      lsf.apply(c,st,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void applyLhs( float[][] wp, float[][] x, float[][] y) {
    zero(y);
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        wpi = (wpi>0.05f)?wpi:0.05f;
        float wps = wpi*wpi;
        float d11 = wps;
        float d22 = wps;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1;
        float y2 = d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  
  //private static void makeRhs
  private static void makeRhs(float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    zero(y);
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        wpi = (wpi>0.05f)?wpi:0.05f;
        float p2i = p2[i2][i1];
        float p3i = p3[i2][i1];
        float b11 = wpi;
        float b22 = wpi;
        float x1 = wpi*p2i;
        float x2 = wpi*p3i;
        float y1 = b11*x1;
        float y2 = b22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
 
  // use 3 points to fit a parabolic curve and find its peak
  private float parabolicPeak(float z, float um, float ui, float up) {
    float a = um-up;
    float b = 2.0f*(um+up)-4.0f*ui;
    return (z+a/b);
  }
  private void checkConstraintForce(float[] k2, float[] k3, float[] cf) {
    int np = k2.length;
    div(cf,max(cf),cf);
    for (int ip=0; ip<np; ++ip) {
      int i2 = (int)k2[ip];
      int i3 = (int)k3[ip];
      System.out.println(" i2="+i2+" i3="+i3+" ad="+cf[ip]);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _weight = 0.5f; // -(fxx+fyy)+weight*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter = 10; // external iterations of surface updating
}
