// SurfaceNets, 10 August 2015, Roberto Toro
// translated from Mikola Lysenko
// http://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	float x,y,z;
} float3D;
typedef struct {
	int a,b,c;
} int3D;
typedef struct {
	int np,nt;
	float3D *p;
	int3D *t;
} Mesh;

int cube_edges[24];
int edge_table[256];
int buffer[4096];

void init_surfacenets(void)
{
	int i,j,p,em,k = 0;
	for(i=0; i<8; ++i) {
		for(j=1; j<=4; j=j<<1) {
			p = i^j;
			if(i <= p) {
				cube_edges[k++] = i;
				cube_edges[k++] = p;
			}
		}
	}
	for(i=0; i<256; ++i) {
		em = 0;
		for(j=0; j<24; j+=2) {
			int a = !(i & (1<<cube_edges[j]));	// was !!, which in js turns into boolean false null, undefined, etc
			int b = !(i & (1<<cube_edges[j+1]));
			em |= a != b ? (1 << (j >> 1)) : 0; // was !==
		}
		edge_table[i] = em;
	}
}
void SurfaceNets(float *data, int *dims, float level, Mesh *mesh, int storeFlag)
{ 
	float3D *vertices=mesh->p;
	int3D *faces=mesh->t;
	int n = 0;
	float x[3];
	int R[3];
	float *grid = (float*)calloc(8,sizeof(float));
	int buf_no = 1;
	int *buffer;
	int buffer_length=0;
	int vertices_length=0;
	int faces_length=0;
	int	i,j,k;

	R[0]=1;
	R[1]=dims[0]+1;
	R[2]=(dims[0]+1)*(dims[1]+1);
 
	if(R[2] * 2 > buffer_length)
		buffer = (int*)calloc(R[2] * 2,sizeof(int));

	for(x[2]=0; x[2]<dims[2]-1; ++x[2])
	{
		n+=dims[0];
		buf_no ^= 1;
		R[2]=-R[2];
		
		int m = 1 + (dims[0]+1) * (1 + buf_no * (dims[1]+1));
		for(x[1]=0; x[1]<dims[1]-1; ++x[1], ++n, m+=2)
		for(x[0]=0; x[0]<dims[0]-1; ++x[0], ++n, ++m)
		{
			int mask = 0, g = 0, idx = n;
			for(k=0; k<2; ++k, idx += dims[0]*(dims[1]-2))
			for(j=0; j<2; ++j, idx += dims[0]-2)      
			for(i=0; i<2; ++i, ++g, ++idx)
			{
				float p = data[idx]-level;
				grid[g] = p;
				mask |= (p < 0) ? (1<<g) : 0;
			}
			if(mask == 0 || mask == 0xff)
				continue;
			int edge_mask = edge_table[mask];
			float3D v = {0.0,0.0,0.0};
			int e_count = 0;
			for(i=0; i<12; ++i)
			{
				if(!(edge_mask & (1<<i)))
					continue;
				++e_count;
				int e0 = cube_edges[ i<<1 ];       //Unpack vertices
				int e1 = cube_edges[(i<<1)+1];
				float g0 = grid[e0];                 //Unpack grid values
				float g1 = grid[e1];
				float t  = g0 - g1;                  //Compute point of intersection
				if(fabs(t) > 1e-6)
					t = g0 / t;
				else
					continue;
				k=1;
				for(j=0; j<3; ++j)
				{
					k=k<<1;
					int a = e0 & k;
					int b = e1 & k;
					if(a != b)
						((float*)&v)[j] += a ? 1.0 - t : t;
					else
						((float*)&v)[j] += a ? 1.0 : 0;
				}
			}
			float s = 1.0 / e_count;
			for(i=0; i<3; ++i)
				((float*)&v)[i] = x[i] + s * ((float*)&v)[i];
			buffer[m] = vertices_length;
			if(storeFlag)
				vertices[vertices_length++]=v;
			else
				vertices_length++;
			for(i=0; i<3; ++i)
			{
				if(!(edge_mask & (1<<i)) )
					continue;
				int iu = (i+1)%3;
				int iv = (i+2)%3;
				if(x[iu] == 0 || x[iv] == 0)
					continue;
				int du = R[iu];
				int dv = R[iv];
				
				if(storeFlag)
				{
					if(mask & 1)
					{
						faces[faces_length++]=(int3D){buffer[m], buffer[m-du-dv], buffer[m-du]};
						faces[faces_length++]=(int3D){buffer[m], buffer[m-dv], buffer[m-du-dv]};
					}
					else
					{
						faces[faces_length++]=(int3D){buffer[m], buffer[m-du-dv], buffer[m-dv]};
						faces[faces_length++]=(int3D){buffer[m], buffer[m-du], buffer[m-du-dv]};
					}
				}
				else
					faces_length+=2;
			}
		}
	}
	mesh->np=vertices_length;
	mesh->nt=faces_length;
}

int main(int argc, char *argv[])
{
//	FILE *f;
	float	*vol;
	int		dim[3];
	int	i,j,k;
	Mesh	m;
	
	dim[0]=100;
	dim[1]=100;
	dim[2]=100;
	vol=(float*)calloc(dim[0]*dim[1]*dim[2],sizeof(float));
	
	for(i=0;i<dim[0];i++)
	for(j=0;j<dim[1];j++)
	for(k=0;k<dim[2];k++)
	if(pow(i-dim[0]/2,2)+pow(j-dim[1]/2,2)+pow(k-dim[2]/2,2)<10*10)
		vol[k*dim[1]*dim[0]+j*dim[0]+i]=100;

	init_surfacenets();
	SurfaceNets(vol,dim,50.0,&m,0);	// 1st pass: evaluate memory requirements
	m.p=(float3D*)calloc(m.np,sizeof(float3D));
	m.t=(int3D*)calloc(m.nt,sizeof(int3D));
	SurfaceNets(vol,dim,50.0,&m,1);	// 2nd pass: store vertices and triangles
	
	printf("%i %i\n",m.np,m.nt);
	for(i=0;i<m.np;i++)
		printf("%f %f %f\n",m.p[i].x,m.p[i].y,m.p[i].z);
	for(i=0;i<m.nt;i++)
		printf("%i %i %i\n",m.t[i].a,m.t[i].b,m.t[i].c);
	
	return 0;

}