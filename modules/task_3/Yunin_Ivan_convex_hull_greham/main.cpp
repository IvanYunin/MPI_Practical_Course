	// Copyright  2018 Ivan Yunin
#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <stack> 
#include <vector> 
using namespace std; 
  
struct Point{
	int x;
	int y;
};
// stores the centre of polygon (It is made 
// global because it is used in compare function) 
Point p0; 
 int f = 0;
void swap(Point &p1, Point &p2) { 
    Point temp = p1; 
    p1 = p2; 
    p2 = temp; 
} 

Point next_to_top(vector <Point> &S) { 
    Point p = S.back(); 
    S.pop_back(); 
    Point res = S.back(); 
    S.push_back(p); 
    return res; 
}   
// Checks whether the line is crossing the polygon 
int orientation(Point a, Point b, 
                Point c) 
{ 
    int res = (b.y-a.y)*(c.x-b.x) - 
              (c.y-b.y)*(b.x-a.x); 
  
    if (res == 0) 
        return 0; 
    if (res > 0) 
        return 1; 
    return -1; 
} 


int distance(Point p1, Point p2) { 
    return (p1.x - p2.x)*(p1.x - p2.x) + 
          (p1.y - p2.y)*(p1.y - p2.y); 
} 

// compare function for sorting 
int compare(const void *vp1, const void *vp2) { 
    Point *p1 = (Point *)vp1; 
    Point *p2 = (Point *)vp2; 
    int o = orientation(p0, *p1, *p2); 
    if (o == 0) 
        if (distance(p0, *p2) >= distance(p0, *p1)) 
	        return -1; 
	    else return 1; 
    if (o == -1)
        return -1;
    else return 1; 
}
  
// Finds upper tangent of two polygons 'a' and 'b' 
// represented as two vectors. //rerurn size of vector
int merger(Point *a, Point *b, Point *out, int size_a, int size_b) { 
    // n1 -> number of points in polygon a 
    // n2 -> number of points in polygon b 
    int n1 = size_a, n2 = size_b; 
  
    int ia = 0, ib = 0; 
    for (int i=1; i<n1; i++) 
        if (a[i].x > a[ia].x) 
            ia = i; 
  
    // ib -> leftmost point of b 
    for (int i=1; i<n2; i++) 
        if (b[i].x < b[ib].x) 
            ib=i; 
  
    // finding the upper tangent 
    int inda = ia, indb = ib; 
    bool done = 0; 
    while (!done) { 
        done = 1; 
        while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >=0) 
            inda = (inda + 1) % n1; 
  
        while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <=0) { 
            indb = (n2+indb-1)%n2; 
            done = 0; 
        } 
    } 
  
    int uppera = inda, upperb = indb; 
    inda = ia, indb=ib; 
    done = 0; 
    int g = 0; 
    while (!done)//finding the lower tangent 
    { 
        done = 1; 
        while (orientation(a[inda], b[indb], b[(indb+1)%n2])>=0) 
            indb=(indb+1)%n2; 
  
        while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1])<=0) 
        { 
            inda=(n1+inda-1)%n1; 
            done=0; 
        } 
    } 
  
    int lowera = inda, lowerb = indb; 
 
    //ret contains the convex hull after merging the two convex hulls 
    //with the points sorted in anti-clockwise order 
    int ind = uppera;
	int out_size=0; 
    out[out_size++]=a[uppera]; 
    while (ind != lowera) 
    { 
        ind = (ind+1)%n1; 
        out[out_size++]=a[ind]; 
    } 
  
    ind = lowerb; 
    out[out_size++]=b[lowerb]; 
    while (ind != upperb) 
    { 
        ind = (ind+1)%n2; 
        out[out_size++]=b[ind]; 
    } 
    return out_size; 
  
} 
  
int convex_hull(Point *points,Point *out, int size){
   int ymin = points[0].y, min = 0; 
   int n = size;	
   for (int i = 1; i < n; i++) { 
     int y = points[i].y; 
     if ((y < ymin) || (ymin == y && 
         points[i].x < points[min].x)) 
        ymin = points[i].y, min = i; 
   } 
   swap(points[0], points[min]);
//   cout << "begin:\n";	
   p0 = points[0]; 
   qsort(&points[1], n-1, sizeof(Point), compare);  
//   cout << "after qsort:\n";		
   int m = 1;
   for (int i=1; i<n; i++) { 
       while (i < n-1 && orientation(p0, points[i], 
                                    points[i+1]) == 0) 
          i++; 
       points[m] = points[i]; 
       m++; 
   } 
//	cout << "afrer eql points:\n";	
   if (m < 3) return 0; 
   int out_size = 0; 
   out[out_size++] =points[0]; 
   out[out_size++] =points[1];
   out[out_size++] =points[2]; 
   for (int i = 3; i < m; i++) { 
      while (orientation(out[out_size-1], out[out_size], points[i]) != -1) 
         out_size--;//S.pop_back(); 
      out[out_size++]=points[i]; 
   } 
//	cout << "after main while:\n";	
   return out_size;
   //while (!S.empty()) { 
   //    Point p = S.back();
   //    cout << "(" << p.x << ", " << p.y <<")" << endl; 
   //    S.pop_back(); 
   //} 
} 

// Returns the convex hull for the given set of points 
//vector<Point> divide( vector<Point> &a) 
//{   
//    // If the number of points is less than 6 then the 
//    // function uses the brute algorithm to find the 
//    // convex hull 
//    if (a.size() <= 50) 
//        return convex_hull(a); 
//  
//    // left contains the left half points 
//    // right contains the right half points 
//    vector<Point>left, right; 
//    for (int i=0; i<a.size()/2; i++) 
//        left.push_back(a[i]); 
//    for (int i=a.size()/2; i<a.size(); i++) 
//        right.push_back(a[i]); 
//    cout << "divide"<<f<<"\n";
//    // convex hull for the left and right sets 
//    vector<Point>left_hull = divide(left); 
//    vector<Point>right_hull = divide(right); 
//    cout << "before merge"<<f<<"\n";
//    // merging the convex hulls 
//    return merger(left_hull, right_hull); 
//} 

int compare_X(const void *vp1, const void *vp2) { 
    Point *p1 = (Point *)vp1; 
    Point *p2 = (Point *)vp2; 
    return ( p1->x - p2->x );
}

void init_map(Point * res,int n,int u_bound,int l_bound){
	Point p;
	for (int i = 0; i < n; i++){
		p.y = 1+rand()%(u_bound-1);
		p.x = 1+rand()%(l_bound-1);
		res[i] = p;
	}
}


int main(int argc, char*argv[]) {
	srand(static_cast<int>(time(0)));
	int num_p = 100;
	int proc_num, proc_id, flag, sub_num_p, seq_res_size, sub_res_size;
	double s_time_start = 0.0, p_time_start = 0.0;
    double s_time_finish = 0.0, p_time_finish = 0.0;
	Point *points, *seq_res, *par_res, *sub_points, *sub_res, *points2;
	Point for_disp;
	MPI_Datatype PNT;
    MPI_Datatype type[2] = { MPI_INT, MPI_INT };
	MPI_Aint disp[2];
	MPI_Aint addresses[3];
    int blocklen[2] = { 1, 1 };
	if (argc > 1) {
        num_p = atoi(argv[1]);
	}
    MPI_Init(&argc, &argv);
    MPI_Initialized(&flag);
    if (!flag) {
        std::cout << "Init MPI Error";
        return 0;
    }
	MPI_Address(&(for_disp), &addresses[0]);
    MPI_Address(&(for_disp.x), &addresses[1]);
    MPI_Address(&(for_disp.y), &addresses[2]);
    disp[0] = addresses[1] - addresses[0];
    disp[1] = addresses[2] - addresses[0];
	
	MPI_Type_create_struct(2, blocklen, disp, type, &PNT);
    MPI_Type_commit(&PNT); 
  	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	//sub_num_p = static_cast<int>(ceil(static_cast<double>(num_p)/
    //                                    (static_cast<double>(proc_num))));
	
	//int hull_sizes[proc_num];
	sub_num_p = num_p/proc_num;
	sub_points = new Point[sub_num_p];
	sub_res = new Point[sub_num_p];
	MPI_Barrier(MPI_COMM_WORLD);
	if (proc_id == 0) {
		points = new Point[num_p];
		points2 = new Point[num_p];
		seq_res = new Point[num_p];
		init_map(points,num_p,num_p,num_p);
		for (int i = 0; i < num_p; i++)
		points2[i] = points[i];
		s_time_start = MPI_Wtime();
    	qsort(points2, num_p, sizeof(Point), compare_X);
  		seq_res_size = convex_hull(points, seq_res, num_p); 
		s_time_finish = MPI_Wtime();
	}
	MPI_Barrier(MPI_COMM_WORLD);
    p_time_start = MPI_Wtime();
	MPI_Scatter(points2 , sub_num_p, PNT, sub_points,
	sub_num_p, PNT, 0, MPI_COMM_WORLD);
	cout<<proc_id<<endl;
	for (int i = 0; i < sub_num_p;  i++) 
      cout << sub_points[i].x << " "
           << sub_points[i].y << endl; 
    cout<<endl;
	sub_res_size = convex_hull(sub_points,sub_res,sub_num_p);
	MPI_Finalize();
   //for (auto e:ans) 
   //   cout << e.x << " "
   //        << e.y << endl; 
	//cout<<endl;
	//for (auto e:ans2) 
      // cout << e.x << " "
        //    << e.y << endl; 
  
    return 0; 
} 