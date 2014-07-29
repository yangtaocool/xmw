package he;

import java.io. *;
import edu.mines.jtk.dsp. *;
import edu.mines.jtk.util. *;
import edu.mines.jtk.interp. *;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

public class SurfaceRefine {
  
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
  // Refines the surface using the active-surface method that honors amplitudes.
  public void surfaceRefine(float[][] scale, float[][][] Eext, float[][] surf, float h){
    int n2 = surf.length; 
    int n1 = surf[0].length; 
    int nz = Eext[0][0].length;
    float[][] hz = new float[2][n2]; float[][] b = new float[n2][n1];
    Sampling sy1 = new Sampling(nz,1.f,0.f); Sampling sx1 = new Sampling(n2,1.,0.); 
    System.out.println("Surface refinment using the active-surface method......");
    for (int n=1; n<_exniter2; n++){
      VecArrayFloat2 vb = new VecArrayFloat2(b);
      VecArrayFloat2 vx = new VecArrayFloat2(surf);
      Smoother2 smoother2 = new Smoother2(_sigma1,_sigma2);
      A a = new A(scale,smoother2,_step,_weight2);
      CgSolver cs = new CgSolver(0.01,_niter);
      vb.zero();
      computeFext(scale,surf,Eext,h,_step,b);
      smoother2.applyTranspose(b);
      cs.solve(a,vb,vx);
      surf = vx.getArray();
      smoother2.apply(surf);
    }
    //return surf;
  }
  
  private static class A implements CgSolver.A {
    A(float[][] scale, Smoother2 s2, float step, float mu){
      _s2 = s2;
      _step = (double)step;
      _mu = (double)mu;
      _scale = scale;
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
      applyLhs(_scale,z,lx);
      applyLhs(_scale,lx,y);
      v2y.add(_step*_mu,v2z,1.0);
      _s2.applyTranspose(y);  
    }
    private double _step;
    private double _mu;
    private Smoother2 _s2;
    private float[][] _scale;
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
  private static void laplacian(float[][] scale,float[][] x, float[][] y){
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
    mul(y,scale,y);
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
  
  
  private static void computeFext(float[][] scale,float[][] surf, float[][][] Eext, float h,float step,float[][] b){
    int n3 = Eext.length;
    int n2 = Eext[0].length; 
    int n1 = Eext[0][0].length; 
    SincInterpolator Si = new SincInterpolator();
    //Si.setUniform(n1,1,0,n2,1,0,n3,1,0,Eext);
    float x,y,z,Ed,Eu,u;
    for (int i3=0; i3<n3; i3++){
      for (int i2=0; i2<n2; i2++){
        x = (float)i3;
        y = (float)i2;
        z = surf[i3][i2];
        Ed = Si.interpolate(n1,1,0,n2,1,0,n3,1,0,Eext,z+h,y,x);
        Eu = Si.interpolate(n1,1,0,n2,1,0,n3,1,0,Eext,z-h,y,x);
        b[i3][i2] = surf[i3][i2]-step*(Ed-Eu)/(2.f*h); 
        //b[i3][i2] = -(Ed-Eu)/(2.f*h); 
      }
    }
    mul(b,scale,b);
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
  private float _weight2 = 50.0f; // weight2*(fxxxx+2fxxyy+fyyyy)+(dA/dz)=0
  private float _mu1 = 1.f;
  private float _mu2 = -0.01f;
  private float _step = 0.5f;
  private float _sigma1 = 0.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 0.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 200; // maximum number of CG iterations
  private int _exniter1 = 0; // external iterations of surface updating
  private int _exniter2 = 80; // external iterations of the active-suface method
  private float _h = 0.05f;//
}
