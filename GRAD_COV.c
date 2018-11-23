
//3 × 3 normal covariance matrix___________________________________________________________________________________


int h = pointhedge(0, @ptnum);


matrix3 test = 0; 
do{
	
   int prev = hedge_prev(0,h);  int pt0 = hedge_srcpoint(0, h);   int pt1 = hedge_dstpoint(0, h);  int pt2 = hedge_srcpoint(0, prev);
         
     vector a = normalize (point( 0 , "P" , hedge_srcpoint(0, h)));
     vector b = normalize (point( 0 , "P" , hedge_dstpoint(0, h))); 
     vector c = normalize (point( 0 , "P" , hedge_srcpoint(0, prev))); 
	
     vector xy = b-a;  vector zy = c-a;
     vector normal = prim( 0 , "normale", hedge_prim( 0, h)); 
	
     float cross = length( cross( xy, zy)); 
     float doat = dot( xy, zy); 
         
     float angle = atan2( cross,doat); 
     3@test += angle * outerproduct( normal , normal); 
         
  h = pointhedgenext(0, h);
	
 }while(h != -1);



//Gradient_______________________________________________________________________________________________
//Appendix B in the Paper 



//__________________________________________________Helper
vector hedge2(int h) {
	vector a = (point(0, "P", hedge_dstpoint(0, h)));
	vector b = (point(0, "P", hedge_srcpoint(0, h)));
	return b - a;
}


void tas(vector normal; float area; int primnum) {
	int points[] = (primpoints(0, primnum));
	int h = primhedge(0, primnum);
	int h2 = hedge_prev(0, h);

	vector a = hedge2(h);
	vector b = hedge2(h2);
	vector g = cross(a, b);
	area = 0.5 * sqrt(dot(g, g));
	normal = normalize(g);

}

vector hedge(int h; int h2) {
	vector a = point(0, "P", h);
	vector b = point(0, "P", h2);
	return a - b;
}
//____________________________________________________FuncDerivative

function vector[] Gradient_InteriorAngle(vector N;   int pt0, pt1, pt2) {
	vector holder[];
	vector vp0, vp1, vp2;
	vector dxh = normalize(hedge(pt0, pt2));
	vector dyh = normalize(hedge(pt1, pt0));
	vector a = cross(N, dxh);
	vector b = cross(N, dyh);
	vp0 = a + b;
	vp1 = b;
	vp2 = a;

	append(holder, vp0);
	append(holder, vp1);
	append(holder, vp2);
	return holder;
}


matrix3 Gradient_Normale(int pt0, pt1; vector normal; float dblarea) {

	vector a = (hedge(pt0, pt1));
	vector croz = cross(a, normal);
	return (outerproduct(croz, normal) / dblarea);
}

//__________________________________________________________________________________Compute_Gradient 
void grad(int ptnum; vector eigenvector) {

	vector grads = set(0, 0, 0);

	int h = pointhedge(0, ptnum);


	float angle = 0;
	float dotpow = 0;
	float dot = 0;
	vector normal;
	float area;

	do {

		int prev = hedge_prev(0, h);
		int pt0 = hedge_srcpoint(0, h);
		int pt1 = hedge_dstpoint(0, h);
		int pt2 = hedge_srcpoint(0, prev);
		int prim = hedge_prim(0, h);

		tas(normal, area, prim);

		float dbl = length(cross(hedge(pt1, pt0), hedge(pt0, pt2)));
		//___________________________EigenVector* Normal   

		dotpow = pow(dot(eigenvector, normal), 2);
		dot = dot(eigenvector, normal);

		//___________________________InteriorAngle(i,j,k)  

		vector ag = normalize(point(0, "P", hedge_srcpoint(0, h)));
		vector bg = normalize(point(0, "P", hedge_dstpoint(0, h)));
		vector cg = normalize(point(0, "P", hedge_srcpoint(0, prev)));
		angle = acos(dot(bg - ag, cg - ag));

		//__________________________GradienentInteriorAngle(i,j,k)  

		vector GradienentInteriorAngle[];
		GradienentInteriorAngle = Gradient_InteriorAngle(
			normal, pt0, pt1, pt2);

		//___________________________GradientNormal(i,j,k) 

		matrix3 pt0mat = Gradient_Normale(pt2, pt1, normal, dbl);
		matrix3 pt1mat = Gradient_Normale(pt0, pt2, normal, dbl);
		matrix3 pt2mat = Gradient_Normale(pt1, pt0, normal, dbl);

		//___________________________Grad(i,j,k); 

		vector a = 0, b = 0, c = 0;
		a += dotpow * GradienentInteriorAngle[0] + 2 * angle * dot * pt0mat * eigenvector;
		b += dotpow * GradienentInteriorAngle[1] + 2 * angle * dot * pt1mat * eigenvector;
		c += dotpow * GradienentInteriorAngle[2] + 2 * angle * dot * pt2mat * eigenvector;


		setpointattrib(0, "gradient2", pt0, a, "set");
		setpointattrib(0, "gradient2", pt1, b, "set");
		setpointattrib(0, "gradient2", pt2, c, "set");
		h = pointhedgenext(0, h);

	} while (h != -1);

}

vector Eigenvector = point(0, "utvector1", @ptnum);
grad(@ptnum, Eigenvector);
