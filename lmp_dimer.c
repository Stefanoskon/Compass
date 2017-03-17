#define cmax_length 1000

#define pi 3.14159265359


#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#include<math.h>

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

int main(int argc, char **argv)
{
	//
	FILE *fp;
	char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length];	
	int i;
	int int_buffer;

	// lammps file variables
	char title[cmax_length];
	int atoms,bonds,angles,dihedrals;
	int atom_types,bond_types,angle_types,dihedral_types;
	int *atom_ID,*molecule_ID,*type,*nx,*ny,*nz;
	int *B_ID,*B_type,*B1,*B2;
	int *A_ID,*A_type,*A1,*A2,*A3;
	int *D_ID,*D_type,*D1,*D2,*D3,*D4;
	int *I_ID,*I_type,*I1,*I2,*I3,*I4;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double *mass;
	double *q,*x,*y,*z;
		
	//
	char command[cmax_length];
	int read_impr;
	int impropers,improper_types;
		
	//
	double xy,xz,yz;
	
	//
	int j;
	int com_atoms,*com_atoms_array;
	double xcom,ycom,zcom;
	double p1x,p1y,p1z;
	double p2x,p2y,p2z;
	double norm2,norm;
	double ux,uy,uz;
	double d;
	double *x_upper,*y_upper,*z_upper;
	double theta;
    
    double shift_x,shift_y,shift_z;
    double shift_x_coord,shift_y_coord,shift_z_coord;
    int shift_atom;
    double shift_theta,shift_d;
	
	//------------------------------------------------------------------
	
	// check if the file contains impropers
	sprintf(command,"grep impropers %s",argv[1]);
	fp=popen(command,"r");
	sprintf(buffer,"%s","\0");	
	fgets(buffer,cmax_length,fp);
	pclose(fp);
	buffer[strcspn(buffer,"\n")]='\0';
	if(strcmp(buffer,"\0")==0){read_impr=0;}else{read_impr=1;}
	
	//------------------------------------------------------------------
	
	getcwd(current_folder,cmax_length);
	
	//
	
	
	sprintf(file_path,"%s/%s",current_folder,"com.dat");
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    com_atoms=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)com_atoms=com_atoms+1;rewind(fp);
    com_atoms_array=(int*)malloc(com_atoms*sizeof(int));
    for(i=0;i<com_atoms;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&com_atoms_array[i]);}
    fclose(fp);

	// read the original lammps data file -- with impropers!!!
	//------------------------------------------------------------------
	
	sprintf(file_path,"%s/%s",current_folder,argv[1]);
	fp=fopen(file_path,"r");
	if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%s",title);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bonds);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angles);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedrals);
	if(read_impr==1){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&impropers);}
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atom_types);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&bond_types);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&angle_types);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&dihedral_types);
	if(read_impr==1){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&improper_types);}
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&xmin,&xmax);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&ymin,&ymax);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf",&zmin,&zmax);
	//fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf",&xy,&xz,&yz);
	
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	
	mass=(double*)malloc((atom_types+1)*sizeof(double));
	for(i=0;i<atom_types;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%lf",&int_buffer,&mass[i]);
	}
	
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	
	atom_ID=(int*)malloc(atoms*sizeof(int));
	molecule_ID=(int*)malloc(atoms*sizeof(int));
	type=(int*)malloc(atoms*sizeof(int));
	nx=(int*)malloc(atoms*sizeof(int));
	ny=(int*)malloc(atoms*sizeof(int));
	nz=(int*)malloc(atoms*sizeof(int));
	q=(double*)malloc(atoms*sizeof(double));
	x=(double*)malloc(atoms*sizeof(double));
	y=(double*)malloc(atoms*sizeof(double));
	z=(double*)malloc(atoms*sizeof(double));
	for(i=0;i<atoms;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d",&atom_ID[i],&molecule_ID[i],&type[i],&q[i],&x[i],&y[i],&z[i],&nx[i],&ny[i],&nz[i]);
	}
	
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	
	B_ID=(int*)malloc(bonds*sizeof(int));
	B_type=(int*)malloc(bonds*sizeof(int));
	B1=(int*)malloc(bonds*sizeof(int));
	B2=(int*)malloc(bonds*sizeof(int));
	for(i=0;i<bonds;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d",&B_ID[i],&B_type[i],&B1[i],&B2[i]);
	}
	
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	
	A_ID=(int*)malloc(angles*sizeof(int));
	A_type=(int*)malloc(angles*sizeof(int));
	A1=(int*)malloc(angles*sizeof(int));
	A2=(int*)malloc(angles*sizeof(int));
	A3=(int*)malloc(angles*sizeof(int));
	for(i=0;i<angles;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d",&A_ID[i],&A_type[i],&A1[i],&A2[i],&A3[i]);
	}
	
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	
	D_ID=(int*)malloc(dihedrals*sizeof(int));
	D_type=(int*)malloc(dihedrals*sizeof(int));
	D1=(int*)malloc(dihedrals*sizeof(int));
	D2=(int*)malloc(dihedrals*sizeof(int));
	D3=(int*)malloc(dihedrals*sizeof(int));
	D4=(int*)malloc(dihedrals*sizeof(int));
	for(i=0;i<dihedrals;++i)
	{
		fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&D_ID[i],&D_type[i],&D1[i],&D2[i],&D3[i],&D4[i]);
	}
	
	if(read_impr==1)
	{
		fgets(buffer,cmax_length,fp);
		fgets(buffer,cmax_length,fp);
		fgets(buffer,cmax_length,fp);
		
		I_ID=(int*)malloc(impropers*sizeof(int));
		I_type=(int*)malloc(impropers*sizeof(int));
		I1=(int*)malloc(impropers*sizeof(int));
		I2=(int*)malloc(impropers*sizeof(int));
		I3=(int*)malloc(impropers*sizeof(int));
		I4=(int*)malloc(impropers*sizeof(int));
		for(i=0;i<impropers;++i)
		{
			fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%d\t%d\t%d",&I_ID[i],&I_type[i],&I1[i],&I2[i],&I3[i],&I4[i]);
		}
	}
	
	fclose(fp);
	
	// add coords unwrap
	
	//
	
	xcom=0.0;ycom=0.0;zcom=0.0;
	for(i=0;i<atoms;++i)
	{
		for(j=0;j<com_atoms;++j)
		{
			if(i+1==com_atoms_array[j])
			{
				xcom=xcom+x[i];ycom=ycom+y[i];zcom=zcom+z[i];
				break;
			}
		}
	}
	xcom=xcom/com_atoms;ycom=ycom/com_atoms;zcom=zcom/com_atoms;
	
	//
	
	
	p1x=x[com_atoms_array[1]-1]-xcom;
	p1y=y[com_atoms_array[1]-1]-ycom;
	p1z=z[com_atoms_array[1]-1]-zcom;
	
	p2x=x[com_atoms_array[0]-1]-xcom;
	p2y=y[com_atoms_array[0]-1]-ycom;
	p2z=z[com_atoms_array[0]-1]-zcom;
	
	//	x	y	z
	//	p1x	p1y	p1z
	//	p2x	p2y	p2z
	
	ux=p1y*p2z-p1z*p2y;
	uy=p1z*p2x-p1x*p2z;
	uz=p1x*p2y-p1y*p2x;
	norm2=ux*ux+uy*uy+uz*uz;
	norm=sqrt(norm2);
	ux=ux/norm;
	uy=uy/norm;
	uz=uz/norm;
	
	d=atof(argv[2]);    // this is the vertical distance
	/*
	printf("%d\n\n",com_atoms*2+2);
	for(i=0;i<com_atoms;++i)printf("C\t%lf\t%lf\t%lf\n",x[com_atoms_array[i]-1],y[com_atoms_array[i]-1],z[com_atoms_array[i]-1]);
	for(i=0;i<com_atoms;++i)printf("C\t%lf\t%lf\t%lf\n",x[com_atoms_array[i]-1]+d*ux,y[com_atoms_array[i]-1]+d*uy,z[com_atoms_array[i]-1]+d*uz);
	printf("point\t%lf\t%lf\t%lf\n",xcom,ycom,zcom);
	printf("point\t%lf\t%lf\t%lf\n",ux+xcom,uy+ycom,uz+zcom);
	*/
	
	x_upper=(double*)malloc(atoms*sizeof(double));
	y_upper=(double*)malloc(atoms*sizeof(double));
	z_upper=(double*)malloc(atoms*sizeof(double));
	
	for(i=0;i<atoms;++i)
	{
		x_upper[i]=x[i]+d*ux;
		y_upper[i]=y[i]+d*uy;
		z_upper[i]=z[i]+d*uz;
		}
	
	xcom=0.0;ycom=0.0;zcom=0.0;
	for(i=0;i<atoms;++i)
	{
		for(j=0;j<com_atoms;++j)
		{
			if(i+1==com_atoms_array[j])
			{
				xcom=xcom+x_upper[i];ycom=ycom+y_upper[i];zcom=zcom+z_upper[i];
				break;
			}
		}
	}
	xcom=xcom/com_atoms;ycom=ycom/com_atoms;zcom=zcom/com_atoms;
    
    //
    
    shift_atom=atoi(argv[4]);   // atom id on the basal molecule to define the shift vector
    shift_theta=atof(argv[5]);  // shift vector rotation angle about ^u
    shift_d=atof(argv[6]);      // shift distance
    
    shift_x_coord=x_upper[shift_atom-1];
    shift_y_coord=y_upper[shift_atom-1];
    shift_z_coord=z_upper[shift_atom-1];
    
    
    
    
    rotate(ux,uy,uz,xcom,ycom,zcom,shift_theta,shift_x_coord,shift_y_coord,shift_z_coord,&shift_x_coord,&shift_y_coord,&shift_z_coord);
    
    
    
    
    
    
    shift_x=shift_x_coord-xcom;
    shift_y=shift_y_coord-ycom;
    shift_z=shift_z_coord-zcom;
    norm2=shift_x*shift_x+shift_y*shift_y+shift_z*shift_z;
    norm=sqrt(norm2);
    shift_x=shift_x/norm;
    shift_y=shift_y/norm;
    shift_z=shift_z/norm;
    
    
    
    
    
    
    
    for(i=0;i<atoms;++i)
    {
        x_upper[i]=x_upper[i]+shift_d*shift_x;
        y_upper[i]=y_upper[i]+shift_d*shift_y;
        z_upper[i]=z_upper[i]+shift_d*shift_z;
    }
    
    
    
    //
	
    xcom=0.0;ycom=0.0;zcom=0.0;
    for(i=0;i<atoms;++i)
    {
        for(j=0;j<com_atoms;++j)
        {
            if(i+1==com_atoms_array[j])
            {
                xcom=xcom+x_upper[i];ycom=ycom+y_upper[i];zcom=zcom+z_upper[i];
                break;
            }
        }
    }
    xcom=xcom/com_atoms;ycom=ycom/com_atoms;zcom=zcom/com_atoms;
    
	theta=atof(argv[3]);    // twist angle
	for(i=0;i<atoms;++i)
		rotate(ux,uy,uz,xcom,ycom,zcom,theta,x_upper[i],y_upper[i],z_upper[i],&x_upper[i],&y_upper[i],&z_upper[i]);
		
	//
		
	printf("%s\n",title);
	printf("\n");
	printf("%d atoms\n",atoms*2);
	printf("%d bonds\n",bonds*2);
	printf("%d angles\n",angles*2);
	printf("%d dihedrals\n",dihedrals*2);
	if(read_impr==1)printf("%d impropers\n",impropers*2);
	printf("\n");
	printf("%d atom types\n",atom_types);
	printf("%d bond types\n",bond_types);
	printf("%d angle types\n",angle_types);
	printf("%d dihedral types\n",dihedral_types);
	if(read_impr==1)printf("%d improper types\n",improper_types);
	printf("\n");
	printf("%lf %lf xlo xhi\n",xmin,xmax);
	printf("%lf %lf ylo yhi\n",ymin,ymax);
	printf("%lf %lf zlo zhi\n",zmin,zmax);
	//printf("%lf %lf %lf xy xz yz\n",xy,xz,yz);
	printf("\n");
	printf("Masses\n");
	printf("\n");
	for(i=0;i<atom_types;++i)
	{
		printf("%d\t%lf\n",i+1,mass[i]);
	}
	printf("\n");
	printf("Atoms\n");
	printf("\n");
	for(i=0;i<atoms;++i)
	{
		printf("%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",atom_ID[i],molecule_ID[i],type[i],q[i],x[i],y[i],z[i],nx[i],ny[i],nz[i]);
	}
	for(i=0;i<atoms;++i)
	{
		printf("%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",atom_ID[i]+atom_ID[atoms-1],molecule_ID[i]+1,type[i],q[i],x_upper[i],y_upper[i],z_upper[i],nx[i],ny[i],nz[i]);
	}
	printf("\n");
	printf("Bonds\n");
	printf("\n");
	for(i=0;i<bonds;++i){printf("%d\t%d\t%d\t%d\n",B_ID[i],B_type[i],B1[i],B2[i]);}
	for(i=0;i<bonds;++i){printf("%d\t%d\t%d\t%d\n",B_ID[i]+B_ID[bonds-1],B_type[i],B1[i]+atoms,B2[i]+atoms);}
	printf("\n");
	printf("Angles\n");
	printf("\n");
	for(i=0;i<angles;++i){printf("%d\t%d\t%d\t%d\t%d\n",A_ID[i],A_type[i],A1[i],A2[i],A3[i]);}
	for(i=0;i<angles;++i){printf("%d\t%d\t%d\t%d\t%d\n",A_ID[i]+A_ID[angles-1],A_type[i],A1[i]+atoms,A2[i]+atoms,A3[i]+atoms);}
	printf("\n");
	printf("Dihedrals\n");
	printf("\n");
	for(i=0;i<dihedrals;++i){printf("%d\t%d\t%d\t%d\t%d\t%d\n",D_ID[i],D_type[i],D1[i],D2[i],D3[i],D4[i]);}
	for(i=0;i<dihedrals;++i){printf("%d\t%d\t%d\t%d\t%d\t%d\n",D_ID[i]+D_ID[dihedrals-1],D_type[i],D1[i]+atoms,D2[i]+atoms,D3[i]+atoms,D4[i]+atoms);}
	if(read_impr==1)
	{
		printf("\n");
		printf("Impropers\n");
		printf("\n");
		for(i=0;i<impropers;++i){printf("%d\t%d\t%d\t%d\t%d\t%d\n",I_ID[i],I_type[i],I1[i],I2[i],I3[i],I4[i]);}
		for(i=0;i<impropers;++i){printf("%d\t%d\t%d\t%d\t%d\t%d\n",I_ID[i]+I_ID[impropers-1],I_type[i],I1[i]+atoms,I2[i]+atoms,I3[i]+atoms,I4[i]+atoms);}
	}
	
	//
	free(atom_ID);free(molecule_ID);free(type);free(nx);free(ny);free(nz);
	free(B_ID);free(B_type);free(B1);free(B2);
	free(A_ID);free(A_type);free(A1);free(A2);free(A3);
	free(D_ID);free(D_type);free(D1);free(D2);free(D3);free(D4);
	
	free(mass);
	free(q);free(x);free(y);free(z);
	
	//
	
	free(com_atoms_array);
	free(x_upper);
	free(y_upper);
	free(z_upper);
}

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new)
{
    double norm2,norm;
    double rads;
    rads=theta*pi/180.0;
    norm2=ux*ux+uy*uy+uz*uz;norm=sqrt(norm2);
    ux=ux/norm;uy=uy/norm;uz=uz/norm;
    *x_new=(ux*ux*(1.0-cos(rads))+cos(rads))*(x_old-x_center)+(ux*uy*(1.0-cos(rads))-uz*sin(rads))*(y_old-y_center)+(ux*uz*(1.0-cos(rads))+uy*sin(rads))*(z_old-z_center)+x_center;
    *y_new=(ux*uy*(1.0-cos(rads))+uz*sin(rads))*(x_old-x_center)+(uy*uy*(1.0-cos(rads))+cos(rads))*(y_old-y_center)+(uy*uz*(1.0-cos(rads))-ux*sin(rads))*(z_old-z_center)+y_center;
    *z_new=(ux*uz*(1.0-cos(rads))-uy*sin(rads))*(x_old-x_center)+(uy*uz*(1.0-cos(rads))+ux*sin(rads))*(y_old-y_center)+(uz*uz*(1.0-cos(rads))+cos(rads))*(z_old-z_center)+z_center;

}

