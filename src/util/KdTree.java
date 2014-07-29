package util;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

// Xinming Wu
// 2013.11.24
public class KdTree {
  /**
   * @param x array[k][n] of coordinates of scattered points.
   * k dimensions
   * n number of scattered points
   * array x is referenced, not copied.
   */
  public KdTree(float[][] x) {
    _x = x;
    _k = x.length;
    _n = x[0].length;
    _i = new int[_n];
    for (int i=0; i<_n; ++i)
      _i[i] = i;
    _ip = 0; //for ploting, label of partition
    _pid = new float[3][_n]; //for ploting, partition ID
    _root = new Node(0,_n-1);//build a k-d tree
  } 
  public void setNLeaf(int nLeaf) {
    _nLeaf = nLeaf;
  }
  // For test only
  public int findNearestSlow(float[] xq) {
    int ni = 0;
    float[] xi = new float[_k];
    for(int j=0; j<_k; ++j) 
      xi[j] = _x[j][ni]; 
    float dmin = distance(xi,xq);
    for (int i=1; i<_n; ++i) {
      for(int j=0; j<_k; ++j) 
        xi[j] = _x[j][i]; 
      float di = distance(xi,xq);
      if(di<dmin) {
        ni = i;
	dmin = di;
      }
    }
    return ni;
  }

  // For test only
  public int[] findInRangeSlow(float[] xmin, float[] xmax) {
    ArrayList<int[]> list = new ArrayList<int[]>();
    for (int i=0; i<_n; ++i) {   
      boolean inside  = true;
      for (int j=0; j<_k && inside; ++j) {   
        float xji = _x[j][i];
	if(xji<xmin[j] || xji>xmax[j])
	  inside = false;
      }
      if(inside) {
        int[] ii = {i};
	list.add(ii);
      }
    }
    return getArray(list);
  }

  public int findNearest(float[] xq) {
    Search search = new Search(xq);
    findNearest(_root,search);
    return search.i;
  }

  public int[] findInRange(float[] xmin, float[] xmax) {
    BoxSearch search = new BoxSearch(xmin,xmax);
    findInRange(_root,search);
    return getArray(search.list);
  }

  /////////////////////////////////////////////////////////////////
  // private
  private int _k;
  private int _n;
  private int[] _i;
  private Node _root;
  private float[][] _x;// the samples, must not change
  private int _nLeaf = 10;
  // for ploting
  public int _ip;
  public float[][] _pid;
  
  private boolean findNearest(Node node, Search search) {
    // if a leaf, check out the points inside a bucket
    if(node.isLeaf) {
      int p = node.p;
      int q = node.q;
      float[] xq = search.x;
      float[] xi = new float[_k];
      for (int i=p; i<=q; ++i) {
	int ii = _i[i];
        for (int j=0; j<_k; ++j)	
	  xi[j] = _x[j][ii];
	float ds = distance(xq,xi);
	if(ds<search.ds) {
          search.i  = ii;
	  search.ds = ds;
	}
      }
      return ballWithinBounds(search);
    }

    int jp = node.jv;
    float vp = node.xv;
    float[] xq   = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    // left son is closer
    if (xq[jp]<=vp) {
      float temp = xmax[jp];
      xmax[jp] = vp;
      if (findNearest(node.left,search)) return true;
      xmax[jp] = temp; 
    } else { // right son is closer
      float temp = xmin[jp]; 
      xmin[jp] = vp;
      if (findNearest(node.right,search)) return true;
      xmin[jp] = temp;
    }

    if (xq[jp]<=vp) {
      float temp = xmin[jp];
      xmin[jp] = vp;
      if (boundsOverlapBall(search) && findNearest(node.right,search))
	return true;
      xmin[jp] = temp;
    } else {
      float temp = xmax[jp];
      xmax[jp] = vp;
      if (boundsOverlapBall(search) && findNearest(node.left,search))
        return true;
      xmax[jp] = temp;
    }
    return ballWithinBounds(search);
  }

  private void findInRange(Node node, BoxSearch search) {
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    // if leaf, check out points in the bucket
    if(node.isLeaf) {
      int p = node.p;
      int q = node.q;
      for (int i=p; i<=q; ++i) {
	boolean inside = true;
        for (int j=0; j<_k && inside; ++j) {
	  float xji = _x[j][_i[i]]; 
          if(xji<xmin[j] || xji>xmax[j]) 
	    inside = false;
	}
	if(inside) {
	  int[] ii = {_i[i]};
	  search.list.add(ii);
	}
      }
    }else {//if not leaf, check out sons 
      int jv = node.jv;
      float xv = node.xv;
      if (xmin[jv]<=xv) findInRange(node.left,search);
      if (xmax[jv]>=xv) findInRange(node.right,search);
    }
  }
  private class Search {
    private int i; // index of closest point
    private float ds; // distance-squared to closest point
    private float[] x;
    private float[] xmin; // current lower bounds during a search
    private float[] xmax; // current upper bounds during a search
    private Search(float[] x) {
      this.i = -1;
      this.x = x;
      ds = Float.MAX_VALUE;
      xmin = new float[_k];
      xmax = new float[_k];
      for (int j=0; j<_k; ++j) {
	xmin[j] = -Float.MAX_VALUE;
	xmax[j] =  Float.MAX_VALUE;
      }
    }
  }

  private class BoxSearch {
    private float[] xmin; // current lower bounds during a search
    private float[] xmax; // current upper bounds during a search
    private ArrayList<int[]> list; // list[n][1]: list of indices for n samples 
    private BoxSearch(float[] xmin, float[] xmax) {
      this.xmin = xmin;
      this.xmax = xmax;
      this.list = new ArrayList<int[]>();
    }
  }

  
  private boolean ballWithinBounds(Search search) {
    float ds = search.ds;
    float[] x = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    for (int j=0; j<_k; ++j) {
      float sd1 = coordinateDistance(j,x,xmin);
      if (sd1<ds) return false;
      float sd2 = coordinateDistance(j,x,xmax);
      if (sd2<ds) return false;
    }
    return true;
  }

  private boolean boundsOverlapBall(Search search) {
    float ds = search.ds;
    float[] x = search.x;
    float[] xmin = search.xmin;
    float[] xmax = search.xmax;
    float sum = 0.0f;
    for (int j=0; j<_k; ++j) {
      float xj = x[j];
      if (xj<xmin[j]) {
        sum += coordinateDistance(j,x,xmin);  
	if (sum>ds) return false;
      }
      if (xj>xmax[j]) {
        sum += coordinateDistance(j,x,xmax);  
	if (sum>ds) return false;
      }
    }
    return true;
  }

  
  private static float distance(float[] x1, float[] x2) {
    float ds = 0.0f; 
    int k = x1.length;
    for (int i=0; i<k; ++i) {
      float di = x1[i]-x2[i];
      ds += di*di;
    }
    return ds;
  }

  private static float coordinateDistance(
    int j, float[] x1, float[] x2) {
      float d = x1[j]-x2[j];
      return d*d;
  }

  private class Node {
    int jv;   // dimension split by this node, if not a leaf
    float xv; // split value, if not a leaf 
    int p,q;
    boolean isLeaf;
    Node left, right;
    Node(int p, int q) {
      if (q-p<_nLeaf) {
        this.p = p;
        this.q = q;
        left = null;
        right = null;
        isLeaf = true;
      } 
      else {
        isLeaf = false;
        int m = (p+q)/2;
        jv = dimensionToDivide(p,q);
        medianSplit(p,q,_x[jv],m);
        xv = _x[jv][_i[m]];
	_pid[0][_ip] = _x[0][_i[m]]; // for plotting
	_pid[1][_ip] = _x[1][_i[m]]; // for plotting
	_pid[2][_ip] = jv; // for plotting
	_ip++;
        left  = new Node(p,m);
        right = new Node(m+1,q);
      }
    }

    // find the widest spread dimension to be divided
    private int dimensionToDivide(int p,int q) {
      int jsmax = 0;
      float xsmax = 0.0f;
      for (int j=0; j<_k; ++j) {
        int ii = p;
        int ix = _i[ii];
        float xminj = _x[j][ix];
        float xmaxj = xminj;
        for (ii=p+1; ii<=q; ++ii) {
	  ix = _i[ii];
          float xtemp = _x[j][ix];
          if (xminj>xtemp) xminj = xtemp;
          if (xmaxj<xtemp) xmaxj = xtemp;
        }
        float xsj = xmaxj-xminj;
        if (xsj>xsmax) {
           jsmax = j;
           xsmax = xsj;
        }
      }
      return jsmax; 
    }
    
    private void medianSplit(int p, int q, float[] x, int m) {
      while (p<q) {
        int ib = p;
        int ie = q;
        float vm = x[_i[m]];
        do {
          while (x[_i[ib]]<vm) ib++;
          while (x[_i[ie]]>vm) ie--;
          if (ib<=ie) swap(ib++,ie--,_i);
        } while (ib<=ie);
        if (ie<m) p=ib;
        if (ib>m) q=ie;
      }
    }

    private void swap(int i1, int i2, int[] x) {
      int xi1 = x[i1];
      x[i1] = x[i2];
      x[i2] = xi1;
    }
  }

  private static int[] getArray(ArrayList<int[]> list) {
    int n = list.size();
    int[] ia = new int[n];
    for (int i=0; i<n; ++i)
      ia[i] = list.get(i)[0];
    return ia;
  }
}

