package he;

import java.io. *;
import edu.mines.jtk.dsp. *;
import edu.mines.jtk.util. *;
import edu.mines.jtk.interp. *;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

public class SurfaceExtractorI {
  
  public void setWeights(double w1, double w2){
    _weight1 = (float)w1;
    _weight2 = (float)w2;
  }
  
  public void setSmoothings(double s1, double s2){
    _sigma1 = (float)s1;
    _sigma2 = (float)s2;
  }
  
  public void setCG(double small, int niter){
    _small = (float) small;
    _niter = niter;
  }
  
  public void setExternalIterations(int niter1, int niter2){
    _exniter1 = niter1;
    _exniter2 = niter2;
  }
  // Initializes a surface using the control points
  public float[][] surfaceInitialization(int n1, int n2, float lmt, float[][] cp){
    float[][] surf = new float[n2][n1]; 
    int n = cp.length; 
    ThinPlateSpline2 tps = new ThinPlateSpline2(n1,n2,cp);
    if (n==1){surf=tps.horizontalSurface();}
    else {
	  ThinPlateSpline2.Coefficients cfs = tps.coefficientsFromPoints();
	  surf = cfs.interpolation();
      }
    for (int i2=0; i2<n2; i2++){
      for (int i1=0; i1<n1; i1++){
        if (surf[i2][i1]<0.f){surf[i2][i1]=0.f;}
        if (surf[i2][i1]>lmt){surf[i2][i1]=lmt;}
      }
    }
    return surf;
  }
  
  // Updates the surface using the seismic normal vectors and control points.
  public float[][] surfaceUpdateFromSlopes
    (float[][][] ep, float[][][] p, float[][][] q,
     float[][] surf, float[][] cp, float lmt, float[][] Ipx, int yt)
  {	
    int n1 = surf[0].length; int n2 = surf.length; 
    int nz = p[0][0].length;
    float[][] b = new float[n2][n1]; float[][] h = new float[2 ][n2];
    Sampling sy1 = new Sampling(nz,1.,0.); Sampling sx1 = new Sampling(n2,1.,0.); 
    for (int k=0; k<n2; k++){h[0][k]=(float)k; h[1][k] = surf[k][yt];}
    //plotFrame(Ipx,h,sy1,sx1,"Iteration "+0,"amplitude");
    System.out.println("Surface updating using seismic normal vectors and control points:");
    for (int n=1; n<=_exniter1; n++){
      System.out.println(" Iteration "+n+"......");
      float[][] p2 = new float[n2][n1]; 
      float[][] p3 = new float[n2][n1]; 
      float[][] wp = new float[n2][n1]; 
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vsurf = new VecArrayFloat2(surf);
      Smoother2 smoother2 = new Smoother2(_sigma1, _sigma2);
      A1 a1 = new A1(smoother2,wp);
      CgSolver cs = new CgSolver(_small, _niter);
      vb.zero();
      findSlopes(p,q,ep,surf,p2,p3,wp);
      makeRhs(wp,p2,p3,b);
      smoother2.applyTranspose(b);
      cs.solve(a1, vb, vsurf);
      surf = vsurf.getArray();
      smoother2.apply(surf);
      surfaceCorrection(lmt,surf,cp);
    }
    return surf;
  }
  // Refines the surface using the active-surface method that honors amplitudes.
  public float[][] surfaceRefine(float[][][] Eext, float[][] surf, float h) {
    int n1 = surf[0].length; int n2 = surf.length; int nz = Eext[0][0].length;
    float[][] hz = new float[2][n2]; float[][] b = new float[n2][n1];
    Sampling sy1 = new Sampling(nz,1.f,0.f); Sampling sx1 = new Sampling(n2,1.,0.); 
    System.out.println("Surface refinment using the active-surface method......");
    for (int n=1; n<_exniter2; n++){
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vx = new VecArrayFloat2(surf);
      Smoother2 smoother2 = new Smoother2(_sigma1,_sigma2);
      A2 a2 = new A2(smoother2,_step,_weight2);
      CgSolver cs = new CgSolver(0.01,_niter);
      vb.zero();
      computeFext(surf,Eext,h,_step,b);
      smoother2.applyTranspose(b);
      cs.solve(a2,vb,vx);
      surf = vx.getArray();
      smoother2.apply(surf);
    }
    return surf;
  }
  
  // Active-surface method
  private float[][] activeSurface(float[][][] Eext, float[][] surf, float h){
    int n1 = surf[0].length; int n2 = surf.length;
    float[][] b = new float[n2][n1];
    for (int n=1; n<2; n++){
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vx = new VecArrayFloat2(surf);
      Smoother2 smoother2 = new Smoother2(_sigma1,_sigma2);
      A2 a2 = new A2(smoother2,_step,_weight2);
      CgSolver cs = new CgSolver(0.001,_niter);
      vb.zero();
      computeFext(surf,Eext,h,_step,b);
      smoother2.applyTranspose(b);
      cs.solve(a2,vb,vx);
      surf = vx.getArray();
      smoother2.apply(surf);
    }
    return surf;
  } 
  private static class A1 implements CgSolver.A{
    A1(Smoother2 s2, float[][] wp){
      _s2=s2;
      _wp=wp;
    }  
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      VecArrayFloat2 v2z = v2x.clone();
      float[][] z = v2z.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] yy = new float[n2][n1];
      float[][] yt = new float[n2][n1];
      VecArrayFloat2 v2yy = new VecArrayFloat2(yy);
      v2y.zero();
      v2yy.zero();
      _s2.apply(z);
      applyLhs(_wp,z,y);
      applyLhs(_wp,z,yt);
      applyLhs(_wp,yt,yy);
      v2y.add(1.f,v2yy,0.f);//-(fxx+fyy)+mu*(fxxxx+fyyyy+2fxxyy)
      _s2.applyTranspose(y);
    }
    private Smoother2 _s2;
    private float[][] _wp;
  }
  
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float step, float mu){
      _s2 = s2;
      _step = (double)step;
      _mu = (double)mu;
    }
    public void apply(Vec vx, Vec vy){
      VecArrayFloat2 v2x = (VecArrayFloat2) vx;
      VecArrayFloat2 v2y = (VecArrayFloat2) vy;
      VecArrayFloat2 v2z = v2x.clone();
      v2y.zero();
      float[][] z = v2z.getArray();
      float[][] y = v2y.getArray();
      int n1 = y[0].length; int n2 = y.length;
      float[][] lx = zerofloat(n1,n2);
      _s2.apply(z);
      laplacian(z,lx);
      laplacian(lx,y);
      v2y.add(_step*_mu,v2z,1.0);
      _s2.applyTranspose(y);  
    }
    private double _step;
    private double _mu;
    private Smoother2 _s2;
  } 
  // Smoother used as a preconditioner 
  private static class Smoother2{
	public Smoother2(float sigma1, float sigma2){
	  _sigma1=sigma1;
	  _sigma2=sigma2;
	}
	public void apply(float[][] x){
	  smooth2(_sigma2,x);
	  smooth1(_sigma1,x);
	}
	public void applyTranspose(float[][] x){
	  smooth1(_sigma1,x);
	  smooth2(_sigma2,x);
	}
    private float _sigma1, _sigma2;  
  }
 
  // Smoothing for dimension 1
  private static void smooth1(float sigma, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n1);
    float[] yt = zerofloat(n1);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }
  // Smoothing for dimension 2
  private static void smooth2(float sigma, float[][] x){
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      lsf.apply(c,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }
  
  // Discrete negative Laplacian operator with neumman boundary conditions.
  private static void laplacian(float[][] x, float[][] y){
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=0; i2<n2; ++i2) {
      int m2 = (i2>0)?i2-1:0;
      for (int i1=0; i1<n1; ++i1) {
        int m1 = (i1>0)?i1-1:0;
        float x1 = x[i2][i1]-x[i2][m1];
        float x2 = x[i2][i1]-x[m2][i1];
        y[i2][i1] += x1+x2;
        y[i2][m1] -= x1;
        y[m2][i1] -= x2;
      }
    }
  }

  private static void applyLhs(
    float[][] wp, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
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

  
  // Discrete nagative Laplacian operator with Derelict boundary conditions
  //private static void laplacian1(float[][] x, float[][] y){

  //}
  //private static void makeRhs
  private static void makeRhs(float[][] wp, float[][] p2, float[][] p3, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
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
  
  
  private void surfaceCorrection(float lmt, float[][] surf, float[][] cp){
    int n1 = surf[0].length; int n2 = surf.length; int n = cp.length;
    float[][] cor = new float[n2][n1];
    float[][] dcp = arrayClone(cp);
    for (int i=0; i<n; i++){
      int i1c = (int) dcp[i][0]; int i2c = (int) dcp[i][1];
      dcp[i][2] -= surf[i2c][i1c];
    }
    ThinPlateSpline2 tps = new ThinPlateSpline2(n1,n2,dcp);
    if (n==1){cor=tps.horizontalSurface();}
    else{
      ThinPlateSpline2.Coefficients cfs = tps.coefficientsFromPoints();
      cor = cfs.interpolation();
    }
    for (int i2=0; i2<n2; i2++){
      for (int i1=0; i1<n1; i1++){
	surf[i2][i1] += cor[i2][i1];
	if (surf[i2][i1]<0.0f){surf[i2][i1]=0.0f;}
	if (surf[i2][i1]>lmt){surf[i2][i1]=lmt;}
      }
    }
  }

  
  private static void findSlopes(float[][][] p, float[][][] q, float[][][] ep, float[][] surf, float[][] p2, float[][] p3,float[][] wp){
    int n2 = surf.length; int n1 = surf[0].length;
    for (int i2=0; i2<n2; i2++){
	  for (int i1=0; i1<n1; i1++){
	    int z = Math.round(surf[i2][i1]);
	    p2[i2][i1] = p[i2][i1][z];
	    p3[i2][i1] = q[i2][i1][z];
	    wp[i2][i1] = ep[i2][i1][z];
	  }
	}
  }
  
  private static void computeFext(
    float[][] surf, float[][][] Eext, float h,float step,float[][] b)
  {
    int n1 = Eext[0][0].length; int n2 = Eext[0].length; int n3 = Eext.length;
    SincInterpolator Si = new SincInterpolator();
    //Si.setUniform(n1,1,0,n2,1,0,n3,1,0,Eext);
    float x,y,z,Ed,Eu,u;
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        x = (float)i3;
        y = (float)i2;
        z = surf[i3][i2];
        //u  = Si.interpolate(z,y,x);
        Ed = Si.interpolate(n1,1,0,n2,1,0,n3,1,0,Eext,z+h,y,x);
        Eu = Si.interpolate(n1,1,0,n2,1,0,n3,1,0,Eext,z-h,y,x);
        b[i3][i2] = surf[i3][i2]-step*(Ed-Eu)/(2.f*h); 
        //b[i3][i2] = -(Ed-Eu)/(2.f*h); 
      }
    }
  }
  
  public static float[][] arrayClone(float[][] a){
    int ny = a.length; int nx = a[0].length;
    float[][] b = new float[ny][nx];
    for (int i=0; i<ny; i++){
	  for (int j=0; j<nx; j++){
	    b[i][j] = a[i][j];
	  }
	}
	return b;
  }	
  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _weight1 = 0.5f; // (fxx+fyy)-weight1*(fxxxx+2fxxyy+fyyyy)=px+qy
  private float _weight2 = 0.1f; // (fxxxx+2fxxyy+fyyyy)+weight2*(dA/dz)=0
  private float _mu1 = 1.f;
  private float _mu2 = -0.01f;
  private float _step = 0.1f;
  private float _sigma1 = 6.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 6.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter1 = 10; // external iterations of surface updating
  private int _exniter2 = 12; // external iterations of the active-suface method
  private float _h = 0.05f;//
}
