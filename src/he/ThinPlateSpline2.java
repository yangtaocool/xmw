package he;

import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;

 /**
  * 2D thin-plate spline interpolation
  * interpolate a surface with minimum curvature from scattered points
  * @author Xinming Wu
  * @version 2012.12.08
  */ 

public class ThinPlateSpline2 {
  public ThinPlateSpline2(int n1, int n2, float[][] cp){
    _n1 = n1;
    _n2 = n2;
    _n = cp.length;
    _cp = cp;//cp[i][0]:xi; cp[i][1]:yi; cp[i][2]:zi
  }
  
  // If there is only one control point, 
  // gets a horizontal surface passing through the control point.
  public float[][] horizontalSurface() {
    float[][] surf = new float[_n2][_n1];
    for (int i2=0; i2<_n2; i2++)
      for (int i1=0; i1<_n1; i1++)
	surf[i2][i1] = _cp[0][2];
    return surf;
  }
  
  public class Coefficients{
    public double[][] coef; //coef=[w1,w2,w3,...,wn,a1,a2,a3]
    private Coefficients(double[][] coef) {
      this.coef = coef;
    }
    // Use the coefficents to interpolate a surface.
    public float[][] interpolation() {
      return apply(coef);
    }
    // Inpterpolates a surface using the interpolation-function coefficients.
    private float[][] apply(double[][] coef) {
      float[][] surf = new float[_n2][_n1];
      double vi, ci1, ci2, ddi, temp;
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<_n1; i1++) {
          vi = coef[_n][0] + coef[_n+1][0]*i1 + coef[_n+2][0]*i2;
	  for (int i=0; i<_n; i++) {
            ci1 = _cp[i][0]; 
            ci2 = _cp[i][1];
	    ddi = (i1-ci1)*(i1-ci1) + (i2-ci2)*(i2-ci2);
	    if (ddi==0.0) temp=0.0;
	    else temp = ddi*Math.log(Math.sqrt(ddi));
            vi += coef[i][0]*temp;
            surf[i2][i1] = (float)vi;
          }
        }
      }
      return surf;
    }
  }
  
  // Get the interpolation coeffients by solving the linear equation Mx=b.
  public Coefficients coefficientsFromPoints(){
    double[][] x = new double[_n+3][1];
    double[][] b = makeRhsVector();
    double[][] M = makeLhsMatrix();
    DMatrix DM = new DMatrix(M);
    DMatrix Db = new DMatrix(b);
    DMatrix Dx = new DMatrix(x);
    DMatrixLud dlu = new DMatrixLud(DM);
    Dx = dlu.solve(Db);
    x = Dx.get();
    return new Coefficients(x);
  }
  
  // Gets the matrix M for Mx=b to compute the coefficients x.
  private double[][] makeLhsMatrix(){
    double[][] M = new double[_n+3][_n+3];
    // Constructs K in M.
    for (int i=0; i<_n; i++){
      for (int j=0; j<_n; j++){
        if (i==j) M[i][j]=0.0;
        else {
          double p1 = pow((_cp[i][0]-_cp[j][0]), 2.0);
          double p2 = pow((_cp[i][1]-_cp[j][1]), 2.0);
          double sp = p1+p2;
          M[i][j] = sp*Math.log(Math.sqrt(sp));
	}
      }
    }
    // Constructs P and P^T in M.
    for (int i=0; i<_n; i++){
      for (int j=_n; j<_n+3; j++){
        if (j==_n){M[i][j]=1.f;M[j][i]=1.f;}
	else if (j==_n+1){
          M[i][j] = (double)_cp[i][0]; 
          M[j][i] = (double)_cp[i][0];
	}
	else{
          M[i][j] = (double)_cp[i][1]; 
          M[j][i] = (double)_cp[i][1]; 
	}
      }
    }
    // Sets zeros in the east-south corner with 3*3 elements in M.
    for (int i=_n; i<_n+3; i++){
      for (int j=_n; j<_n+3; j++){
        M[i][j] = 0.0;
      }
    }
    return M; 
  }
  // Gets the matrix b for Mx=b to compute the coefficients x.
  private double[][] makeRhsVector(){
    double[][] b = new double[_n+3][1];
    for (int i=0; i<_n; i++){b[i][0]=(double)_cp[i][2];}
    b[_n  ][0]=0.0;
    b[_n+1][0]=0.0;
    b[_n+2][0]=0.0; 
    return b;
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // private
  private int _n1;// the number of points to be interpolated in the 1st dim
  private int _n2;// the number of points to be interpolated in the 2nd dim
  private int _n; // the number of scattered control points
  private float[][] _cp; // positions of the control points
}
