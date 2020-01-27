//__________________________________________________Helper
vector hedge2( int h){
vector a = (point( 0 , "P" , hedge_dstpoint(0, h))); 
vector b = (point( 0 , "P" , hedge_srcpoint(0, h)));
return a-b; 
}


void tas( vector normal ; float area; int primnum) {
int points[] = (primpoints(0,primnum ));
int h =         primhedge(0,primnum);
int h2 =        hedge_prev( 0, h); 

vector a =  hedge2( h); 
vector b =  hedge2( h2); 
vector g= cross( a,b); 
area= 0.5 * sqrt(dot(g,g)); 
normal =normalize(g); 

}



vector hedge( int h;int h2){
vector a = point( 0 , "P" , h); 
vector b = point( 0 , "P" ,h2);
return a-b; 
}

_____________________________________________________________________

function vector[] Gradient_InteriorAngle( vector N;   int pt0 ,pt1 , pt2  ){
vector holder[]; 
vector vp0, vp1 , vp2; 
vector dxh =(hedge(pt0 , pt1)); 
vector dyh =(hedge(pt2 , pt0));
vector a = cross( (N), dxh); 
vector b = cross( (N), dyh); 
vp0 = a+b ;   
vp1 = a;  
vp2 = b;   

append( holder,-vp0); 
append( holder,vp1); 
append( holder,vp2); 
return holder;  
}


matrix3 Gradient_Normale (int pt0 ,pt1 ; vector normal  ;float dblarea) {

           vector a = (hedge(pt0, pt1) );
           vector croz = cross(a, normal); 
           return ( outerproduct(croz,normal )/dblarea);
} 

//__________________________________________________________________________________Compute_Gradient 


vector grads = set(0,0,0);

int h = pointhedge(0, @ptnum);

vector Eigenvector = point( 0 , "utvector1",@ptnum); 

float angle = 0; 
float dotpow = 0;
float dot = 0;
vector normal; 
float area; 
vector b=set(0,0,0);vector a=set(0,0,0);vector c=set(0,0,0);
          
do{
     
         int prev = hedge_prev(0,h); 
         int pt0 =  hedge_srcpoint(0, h); 
         int pt1 =  hedge_dstpoint(0, h);
         int pt2 =  hedge_srcpoint(0, prev);
         int prim = hedge_prim(0, h); 
         
          tas(normal, area,prim);  
         
      
      
       
        
        
//___________________________EigenVector* Normal   

         dotpow= pow(dot( Eigenvector, normal),2); 
         dot= dot( Eigenvector, normal);  
                                     
//___________________________InteriorAngle(i,j,k)  

         vector ag = (point( 0 ,"P" ,hedge_srcpoint(0, h))); 
         vector bg = (point( 0 , "P" , hedge_dstpoint(0, h))); 
         vector cg = (point( 0 , "P" , hedge_srcpoint(0, prev))); 
         vector  xy = bg-ag; 
         vector  zy = cg-ag; 
         float cross = length( cross( xy, zy)); 
         float doat = dot( xy, zy); 
         angle = atan2( cross,doat); 
         float dbl = length( cross); 
         
         
         
         
//__________________________GradienentInteriorAngle(i,j,k)  

         vector GradienentInteriorAngle[]; 
         GradienentInteriorAngle =   Gradient_InteriorAngle( 
                                     normal, pt0, pt1 , pt2);                 
          
            
//___________________________GradientNormal(i,j,k) 

           matrix3 pt0mat =   Gradient_Normale( pt1,pt2, normal,dbl); 
           matrix3 pt1mat =   Gradient_Normale( pt2,pt0, normal,dbl); 
           matrix3 pt2mat =   Gradient_Normale( pt0,pt1, normal,dbl);  
           
//___________________________Grad(i,j,k)
	
            a = dotpow* GradienentInteriorAngle[0]+2 * angle * dot *    Eigenvector * pt0mat; 
            b = dotpow*GradienentInteriorAngle[1]+2 * angle * dot *    Eigenvector * pt1mat; 
            c = dotpow* GradienentInteriorAngle[2]+2 * angle * dot *     Eigenvector * pt2mat; 
     
              
         
setpointattrib( 0 , "gradient" ,pt0, a, "add");  
setpointattrib( 0 , "gradient" ,pt1, b, "add"); 
setpointattrib( 0 , "gradient" ,pt2, c, "add"); 
h = pointhedgenext(0, h);        
              
}while(h != -1);

