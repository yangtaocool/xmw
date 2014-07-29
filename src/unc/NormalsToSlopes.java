/****************************************************************************
Copyright (c) 2007, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package unc;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Local estimates of orientations of features in images.
 * Methods of this class can compute for each image sample numerous
 * parameters related to orientation. All orientation information 
 * is derived from eigenvectors and eigenvalues of the structure tensor
 * (also called the "gradient-squared tensor"). This tensor is equivalent 
 * to a matrix of 2nd partial derivatives of an autocorrelation evaluated 
 * at zero lag. In other words, orientation is here determined by the 
 * (2-D) ellipse or (3-D) ellipsoid that best fits the peak of the 
 * autocorrelation of image samples in a local window.
 * <p>
 * The coordinate system for a 2-D image has two orthogonal axes 1 and 2, 
 * which correspond to the 1st and 2nd indices of the array containing 
 * image samples. For 2-D images, the eigenvectors are the unit vectors 
 * u = (u1,u2) and v = (v1,v2). The 1st eigenvector u is perpendicular 
 * to the best fitting line, and the 1st component u1 of u is always 
 * non-negative. The 2nd eigenvector v is perpendicular to u such that 
 * the cross product u1*v2-u2*v1 = 1; that is, v1 = -u2 and v2 = u1. 
 * The angle theta = asin(u2) is the angle measured counter-clockwise 
 * between the 1st eigenvector u and axis 1; -pi/2 &lt;= theta &lt;= pi/2.
 * <p>
 * The coordinate system for a 3-D image has three orthogonal axes 1, 2 
 * and 3, which correspond to the 1st, 2nd and 3rd indices of the array 
 * containing image samples. For 3-D images, the eigenvectors are unit 
 * vectors u = (u1,u2,u3), v = (v1,v2,v3), and w = (w1,w2,w3). The 1st 
 * eigenvector u is orthogonal to the best fitting plane, and the 1st 
 * component u1 of u is always non-negative. The 2nd eigenvector v is 
 * orthogonal to the best fitting line within the best fitting plane.
 * The 3rd eigenvector w is orthogonal to both u and v and is aligned
 * with the direction in which the images changes least. The dip angle 
 * theta = acos(u1) is the angle between the 1st eigenvector u and axis 1; 
 * 0 &lt;= theta &lt;= pi/2. The azimuthal angle phi = atan2(u3,u2)
 * is well-defined for only non-zero theta; -pi &lt;= phi &lt;= pi.
 * <p>
 * The local linearity or planarity of features is determined by the
 * eigenvalues. For 2-D images with eigenvalues eu and ev (corresponding 
 * to the eigenvectors u and v), linearity is (eu-ev)/eu. For 3-D
 * images with eigenvalues eu, ev, and ew, planarity is (eu-ev)/eu
 * and linearity is (ev-ew)/eu. Both linearity and planarity are
 * in the range [0,1].
 *
 * @author Xinming Wu, Colorado School of Mines
 * @version 2013.12.07
 */
public class NormalsToSlopes {

  public void setLimits(float pmax) {
    setLimits(-pmax,pmax,-pmax,pmax);
  }
  public void setLimits(float p2min, float p2max, float p3min, float p3max) {
    _p2min = p2min;   
    _p2max = p2max;   
    _p3min = p3min;   
    _p3max = p3max;   
  }

  public void apply(
    float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] p2, float[][][] p3) 
  {
    int n3 = u1.length; 
    int n2 = u1[0].length; 
    int n1 = u1[0][0].length; 
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = u1[i3][i2][i1];
          float u2i = u2[i3][i2][i1];
          float u3i = u3[i3][i2][i1];
	        if (-u2i<_p2min*u1i) u2i = -_p2min*u1i;
          if (-u2i>_p2max*u1i) u2i = -_p2max*u1i;
          if (-u3i<_p3min*u1i) u3i = -_p3min*u1i;
          if (-u3i>_p3max*u1i) u3i = -_p3max*u1i;
          if (u1i==0.0f) {
            p2[i3][i2][i1] = (u2i<0.0f)?_p2max:_p2min;
            p3[i3][i2][i1] = (u3i<0.0f)?_p3max:_p3min;
          } else {
            p2[i3][i2][i1] = -u2i/u1i;
            p3[i3][i2][i1] = -u3i/u1i;
          }
	}
      }
    }
  }
  private float _p2min = -100.0f;
  private float _p2max =  100.0f;
  private float _p3min = -100.0f;
  private float _p3max =  100.0f;
}







