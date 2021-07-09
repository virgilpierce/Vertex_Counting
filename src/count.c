/* 
Program for counting ribbon graphs with vertices of uniform valence.  

Based upon the permutation structure giving vertices and edges.  The pertinent fact form which is that the faces are then give by the permutation:  vertex permutation * edge permuatation.  

See my unpublished paper on the algorithm on the arxive for further references:  http://arxiv.org/abs/math/0610586  
*/
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<ostream>
#include<iostream>
#include<iomanip>
#include<time.h>
using namespace std;


#define vertices 5  /* Number of vertices */
#define degree 4   /* degree of each vertex */

#define arrows 20  /* arrows are half-edges; so this number should be
			vertices * degree
*/

/* The labeling and vertex permutation are choosen such that it is 
vertices = (0 1 ... degree-1) (degree degree+1 ... 2*degree -1) ...
*/

/* The data of the ribbon graph is represented as a pair of permutations.  The faces, are then given by the product of this pair.  The vertex permutation is fixed as above; and then we will cycle through all possible edge permutations (same as Wick pairings).

This is the data structure of a permuatation.
*/
typedef struct
{ int entry[arrows]; } permutation;  

/* The program counts the number of ribbon graphs and bins the results by genus in a data structure of this type */
typedef struct
{ long double entry[10]; } genus_count;


/* check(edge_map) checks that the ribbon graph corresponding to edge_map is connected.  Returns a 1 if it is, a zero if it isn't */
int check(permutation edge_map)
{
    int ans = 1;
    int vertex[vertices];
    
    vertex[0] = 1;
    for( int m=1; m<vertices; m++)
	vertex[m] = 0;

/* here a zero means that the edge_map has not yet visited vertex m.
*/

    for(int k = 0; k<vertices-1; k++)
    {
   	ans = 1;
	for(int m1=0; m1<vertices; m1++)
	/* It is a breadth-first algorithm:  go to a vertex we have visited; and mark of each vertex connected to it as visited if we have not been there already */
	    if(vertex[m1] != 0)
	    {	

		for(int k2=degree*m1; k2<degree*m1+degree; k2++)
		{
		    for(int m2=1; m2 < m1; m2++)
			if( edge_map.entry[k2]/degree == m2)
			    vertex[m2] = 1;
		    for(int m2=m1+1; m2<vertices; m2++)
			if( edge_map.entry[k2]/degree == m2)
			    vertex[m2] = 1;
		}
	    }


	
	for(int k1=0; k1<vertices; k1++)
	    if(vertex[k1] == 0)
		ans =0;
	/* change ans to 0 if there is a vertex left to visit. */

    	if(ans==1)
	   k= vertices-1;

	/* if all vertices have been visited leave the loop over k early */ 
       }


    return ans;
}
    

/* face_count(vertex_map, edge_map) returns the number of faces of the ribbon graph.  The trick is to computer the permutation vertex_map * edge_map; and count the number of orbits. 
*/ 
int face_count(permutation vertex_map, permutation edge_map)
{


    int ans = 0;

    permutation travel;

    for(int k=0; k<arrows; k++)
	travel.entry[k] = 0;

	/* travel represents the arrows of the graph.  We will follow each complete loop of arrows leaving a vertex.  When the loop is complete we have encircled (and thus identified) one face.  */

    for(int k=0; k<arrows; k++)
    {
	int l = k;
	while(travel.entry[k] == 0)
	{
	  /* this gives the permuatation product vertex * edge. */
	    l = vertex_map.entry[ edge_map.entry[l] ];

	/* thus we have arrived at arrow l; and so we mark it as having been used with a 1 */ 
	    travel.entry[l] = 1;

	/* when we return to k we have encircled one face; and so we increment our count by 1 */
	    if( l==k)
		ans ++;
	}
	    
    }



    return ans;
}
	    

/* count computes the genus of the ribbon graph given by vertex_map and edge_map and increments the bin genus_count.
*/
void count(permutation vertex_map, permutation edge_map, genus_count * genus)
{
 

    int conn = check(edge_map);

    if(conn != 0)
    {
	int faces = face_count(vertex_map, edge_map); // compute the number of faces
	int chi = (vertices - arrows/2 + faces); // compute the Euler characteristic
	(*genus).entry[(2 - chi)/2]+= 1.0;  // bin the result by genus.
    }
}


/* The algorithm is an improvement over the gap algorithm in that when changing the number of vertices and the degree nothing is changed beyond those parameters.  The reason is that the for loops cycling over all possible Wick pairings in the gap code have been replaced by nested functions which self terminate.  */

void nest2(int nest_count, permutation *which, permutation *chooser, permutation *edge_map, permutation vertex_map, genus_count *genus)
{
 


    if(nest_count < arrows/2-1 )
    {
	(*edge_map).entry[(*chooser).entry[0]] = (*chooser).entry[(*which).entry[nest_count+1]];
	(*edge_map).entry[(*chooser).entry[(*which).entry[nest_count+1]]] = (*chooser).entry[0];

	for(int j=0; j<arrows - 2*nest_count-1; j++)
	{
	    (*chooser).entry[j] = (*chooser).entry[j+1];
	}
	for(int j=(*which).entry[nest_count+1]-1; j<arrows -2*nest_count-2; j++)
	{
	    (*chooser).entry[j] = (*chooser).entry[j+1];
	}
	

	nest2(nest_count+1, which, chooser, edge_map, vertex_map, genus);
    }
    else
    {
	(*edge_map).entry[(*chooser).entry[0]] = (*chooser).entry[1];
	(*edge_map).entry[(*chooser).entry[1]] = (*chooser).entry[0];
	count(vertex_map, *edge_map, genus);
    }
}

void nest(int nest_count, permutation *which, permutation *chooser, permutation *edge_map, permutation vertex_map, genus_count *genus)
{
 

    if( nest_count < arrows/2 -1)
	for(int i=1; i<arrows - 2*nest_count; i++)
	{
	    (*which).entry[nest_count+1] = i;
	    nest(nest_count+1, which, chooser, edge_map, vertex_map, genus);
	}
    else
    {
	for(int j=0; j<arrows; j++)
	    (*chooser).entry[j] = j;
	nest2(0, which, chooser, edge_map, vertex_map, genus);
    }
}
	
/* nest and nest2 are the nested functions which produce the loops over all possible wick pairings and then bin the result.  

nest builds the chooser; 

nest2 then assigns the elements of chooser to permutation edge_map; when it reaches the end it then calls the check and genus_count functions.

There is opportunity for parallization here:  each call of the nesting function should (in an ideal system) spawn a new process.  

*/


int main()
{

    time_t now;
    time_t start = time(&now);

    permutation vertex_map; //vertex data
    permutation edge_map;  // edge data
    permutation chooser;  // chooser to be used in cycling over all Wick pairings
    permutation which;
    genus_count genus;  // results are binned in genus
    int nest_count = 0; // checks whether we have reached the end of the nested loops.


   // initialize the genus count.
    for(int i1 = 0; i1 < 10; i1++)
	genus.entry[i1] = 0.0;

  // define the vertex permutation.
    for(int i1 =0; i1 < vertices; i1++)
    {
	for( int i2=0; i2<degree-1; i2++)
	{
	    vertex_map.entry[degree*i1 + i2] = degree*i1+i2+1;
	}
	vertex_map.entry[degree*i1 +degree-1] = degree*i1;
    }

   // start the nesting loops.
    nest(0, &which, &chooser, &edge_map, vertex_map, &genus);

  // Display the results.
    for(int i1 = 0; i1<10; i1++)
    {
	cout << "Genus:  " << i1 << "\n";
	    cout << "Number counted:  " << setprecision(60) 
		 << genus.entry[i1] << "\n";
    }

  // Report the time used.  
    time_t ends=time(&now);
    cout << "\n Running Time : " 
	<< difftime(ends, start) << endl;

    return 0;
}
