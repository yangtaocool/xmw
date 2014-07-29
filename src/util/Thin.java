/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package util;

import vec. *;
import java.util.logging.Logger;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Thin {
  
  public void applyForHorizontal(float[][][] u, float[][][] ut1, float[][][] ut2) {
    int n3 = u.length;
    int n2 = u[0].length;
    int n1 = u[0][0].length;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(2.0f);
    ref.apply2(u,u);
    ref.apply2(u,u);
    ref.apply3(u,u);
    ref.apply3(u,u);
    for (int i3=0; i3<n3; ++i3) 
      for (int i1=0; i1<n1; ++i1) 
        for (int i2=1; i2<n2-1; ++i2) {
          float ui =   u[i3][i2  ][i1];
          float um =   u[i3][i2-1][i1];
          float up =   u[i3][i2+1][i1];
          float dumi = -um + ui; 
          float dupi = -up + ui; 
          float sum1 = dumi+dupi;
          float sum2 = abs(dumi)+abs(dupi);
          if (sum1==sum2) {
            ut1[i3][i2  ][i1] = ui;
            ut2[i3][i2  ][i1] = ui;
            ut2[i3][i2-1][i1] = um;
            ut2[i3][i2+1][i1] = up;
          }
       }

   for (int i2=0; i2<n2; ++i2) 
      for (int i1=0; i1<n1; ++i1) 
        for (int i3=1; i3<n3-1; ++i3) {
          float uti = ut2[i3][i2][i1];
          float utm = ut2[i3][i2][i1];
          float utp = ut2[i3][i2][i1];
          float ui = u[i3  ][i2][i1];
          float um = u[i3-1][i2][i1];
          float up = u[i3+1][i2][i1];
          float dumi = -um + ui; 
          float dupi = -up + ui; 
          float sum = uti+utm+utp;
          float sum1 = dumi+dupi;
          float sum2 = abs(dumi)+abs(dupi);
          if (sum1==sum2 && uti>0.0f) {
            ut1[i3  ][i2][i1] = ui;
            ut2[i3  ][i2][i1] = ui;
            ut2[i3-1][i2][i1] = um;
            ut2[i3+1][i2][i1] = up;
          }
       }
  }
}
