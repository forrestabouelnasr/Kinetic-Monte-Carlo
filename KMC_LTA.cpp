/*
 *  Created by forrest on 4/12/10.
 *  
 *  This program simulates the movements of methane at high loadings in 
 *  zeolite LTA at 300K. Hop rate values are taken from TST analysis 
 *  of all-atom simulations. The program executes several kinetic Monte
 *  Carlo simulations and reports the self- and collective-diffusion
 *  coefficients. 
 */
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

double D_self_global;
double D_collective_global;
int n=3;
int n_particles_per_cage=9;
double t_increment=1e-11;
double Hoprates[][3]={
	0, 2.13e11, 0,
	2.96e11, 0, 1.34e11,
	0, 2.91e11, 0.652e11};

int KMC (void)
{
	srand((unsigned)time(0));
    
	n=2;
	int n_cages=15*n*n*n;
	double a, b, b2, c, c2;
	a=12.305/2;
	b=3.9469;
	b2=12.305-b;
	c=1.7564;
	c2=12.305-c;
	double template1[][3] = {//this is one LTA cage, made of 15 subcages
		c, a, a,
		b, b2, b2,
		b, b2, b,
		b, b, b2,
		b, b, b,
		a, c2, a,
		a, a, c2,
		a, a, a,
		a, a, c,
		a, c, a,
		b2, b2, b2,
		b2, b2, b,
		b2, b, b2,
		b2, b, b,
		c2, a, a};
	int ConnectionTypes[][3] = {//define connection types between cage types
		0,1,0,
		1, 0, 3,
		0, 3, 4};
	double template1_dimensions[]={12.305, 12.305, 12.305};
	int template1_CageTypes[]={3, 2, 2, 2, 2, 3, 3, 1, 3, 3, 2, 2, 2, 2, 3};
	for (int k=0;k<15;k++){
		template1_CageTypes[k]--;
	}
	double CageCenters[15*n*n*n][3];
	int CageTypes[n*n*n*15];
	int x, y, z, cage;
	for(x=0;x<n;x++){
		for(y=0;y<n;y++){
			for(z=0;z<n;z++){
				for(cage=0;cage<15;cage++){//create cubic copies of template
					CageCenters[cage+15*z+15*n*y+15*n*n*x][0]=x*template1_dimensions[0]+template1[cage][0];
					CageCenters[cage+15*z+15*n*y+15*n*n*x][1]=y*template1_dimensions[1]+template1[cage][1];
					CageCenters[cage+15*z+15*n*y+15*n*n*x][2]=z*template1_dimensions[2]+template1[cage][2];
					CageTypes[cage+15*z+15*n*y+15*n*n*x]=template1_CageTypes[cage];
				}
			}
		}
	}
	int cage_i, cage_j, CageConnectivities[15*n*n*n][15*n*n*n], i,j;
	int CageConnectivities2[15*n*n*n][15*n*n*n][3];
	double j_coordinates[3], i_coordinates[3], distance_squared;
	for(cage_i=0;cage_i<15*n*n*n;cage_i++){
		for(cage_j=0;cage_j<15*n*n*n;cage_j++){
			CageConnectivities[cage_i][cage_j]=0;
			CageConnectivities2[cage_i][cage_j][0]=0.0;
			CageConnectivities2[cage_i][cage_j][1]=0.0;
			CageConnectivities2[cage_i][cage_j][2]=0.0;
		}
	}
	for(cage_i=0;cage_i<15*n*n*n-1;cage_i++){
		for(cage_j=cage_i+1;cage_j<15*n*n*n;cage_j++){
			for(x=-1;x<=2;x++){
				for(y=-1;y<=2;y++){
					for(z=-1;z<=2;z++){
						j_coordinates[0]=CageCenters[cage_j][0]+double(x)*template1_dimensions[0]*double(n);
						j_coordinates[1]=CageCenters[cage_j][1]+double(y)*template1_dimensions[1]*double(n);
						j_coordinates[2]=CageCenters[cage_j][2]+double(z)*template1_dimensions[2]*double(n);
						i_coordinates[0]=CageCenters[cage_i][0];
						i_coordinates[1]=CageCenters[cage_i][1];
						i_coordinates[2]=CageCenters[cage_i][2];
						distance_squared=
						(j_coordinates[0]-i_coordinates[0])*(j_coordinates[0]-i_coordinates[0])+
						(j_coordinates[1]-i_coordinates[1])*(j_coordinates[1]-i_coordinates[1])+
						(j_coordinates[2]-i_coordinates[2])*(j_coordinates[2]-i_coordinates[2]);
						if(distance_squared<=30){//if cages are close enough, define them as connected
							i=CageTypes[cage_i];
							j=CageTypes[cage_j];
							CageConnectivities[cage_i][cage_j]=ConnectionTypes[i][j];
							CageConnectivities[cage_j][cage_i]=ConnectionTypes[i][j];
							CageConnectivities2[cage_i][cage_j][0]=x;//include PBC
							CageConnectivities2[cage_j][cage_i][0]=-x;
							CageConnectivities2[cage_i][cage_j][1]=y;
							CageConnectivities2[cage_j][cage_i][1]=-y;
							CageConnectivities2[cage_i][cage_j][2]=z;
							CageConnectivities2[cage_j][cage_i][2]=-z;
						}
					}
				}
			}
		}
	}
	for(i=0; i<15*n*n*n; i++){//check for errors
		for(j=0; j<15*n*n*n; j++){
			if(CageConnectivities[i][j]!=CageConnectivities[j][i])
				cout<<"Problem1";
			if(CageConnectivities2[i][j][0]!=-CageConnectivities2[j][i][0])
				cout<<"Problem2";
			if(CageConnectivities2[i][j][1]!=-CageConnectivities2[j][i][1])
				cout<<"Problem3";
			if(CageConnectivities2[i][j][2]!=-CageConnectivities2[j][i][2])
				cout<<"Problem4";
			if(CageConnectivities2[i][j][0]!=CageConnectivities2[i][j][0]*(0.0+(CageConnectivities[i][j])>0))
				cout<<"Problem5 ";
			if(CageConnectivities2[i][j][1]!=CageConnectivities2[i][j][1]*(0.0+(CageConnectivities[i][j])>0))
				cout<<"Problem6 ";
			if(CageConnectivities2[i][j][2]!=CageConnectivities2[i][j][2]*(0.0+(CageConnectivities[i][j])>0))
				cout<<"Problem7 ";
		}
	}
	double CageDimensions[3]={template1_dimensions[0]*n, template1_dimensions[1]*n, template1_dimensions[2]*n};
	double BoxDimensions[3]={CageDimensions[0], CageDimensions[1], CageDimensions[2]};
	for(i=0;i<n_cages;i++){
		for(j=0;j<n_cages;j++){
			if(CageConnectivities[i][j]!=CageConnectivities[j][i]){
				cout<<"Problem A. i= "<<i<<" j= "<<j<<endl;
				return 0;
			}
			if(CageConnectivities[i][j]==0){
				if(CageConnectivities2[i][j][0]!=0){
					cout<<"Problem B 0. i= """<<i<<" j= "<<j<<endl;
					return 0;
				}
				if(CageConnectivities2[i][j][1]!=0){
					cout<<"Problem B 1. i= """<<i<<" j= "<<j<<endl;
					return 0;
				}	
				if(CageConnectivities2[i][j][2]!=0){
					cout<<"Problem B 2. i= "<<i<<" j= "<<j<<endl;
					return 0;
				}	
			}
			if(CageConnectivities[i][j]!=0){
				if(Hoprates[CageTypes[i]][CageTypes[j]]==0){
					cout<<"Problem C. i= "<<i<<" j= "<<j;
					cout<<" CageTypes[i]= "<<CageTypes[i]<<" cagetypes[j]= "<<CageTypes[j]<<endl;
				}
			}
		}
	}
	cout << "Lattice construction complete" << endl;

	int n_particles=n_particles_per_cage*n*n*n;
	double t_final=t_increment*100;
	double t_initialize=t_increment*10;
	//--------------------begin actual simulation----------------------
	double t_record[int(t_final/t_increment)], t=-t_initialize;
	t_record[0]=0;
	int done, selected_particle, selected_cage;
	int CageOccupancy[n_cages];
	for(i=0;i<n_cages;i++){
		CageOccupancy[i]=0;
	}
	double MoleculeDescription[n_particles][6];
	int MoleculeDescription2[n_particles];
	//place particles randomly into cages
	if(n_particles>n_cages){
		cout << "too many particles!"<<endl;
		return 0;
	}
	for(i=0;i<n_particles;i++){
		done=0;
		while(done==0){
			cage=rand()%n_cages;
			if(CageOccupancy[cage]==0){
				done=1;
				MoleculeDescription[i][0]=CageCenters[cage][0];
				MoleculeDescription[i][1]=CageCenters[cage][1];
				MoleculeDescription[i][2]=CageCenters[cage][2];
				MoleculeDescription[i][3]=0;
				MoleculeDescription[i][4]=0;
				MoleculeDescription[i][5]=0;
				MoleculeDescription2[i]=cage;
				CageOccupancy[cage]++;
			}
		}
	}
	double Hopping_Rate_Now[n_particles][n_cages];
	double Hopping_Rate_Cumulative[n_particles][n_cages];
	double Hopping_Rate_Total, tau;
	double t_record_end;
	int particle, r,v;
	double xyz_record[n_particles][7][int(t_final/t_increment)];
	double t_record_max=0;
	for(particle=0;particle<n_particles;particle++){
		for(r=0;r<7;r++){
			xyz_record[particle][r][0]=0;
		}
	}
	int upper_bound, lower_bound, target, record_index=0;
	int p;
	while(t<t_final){
		for(i=0;i<n_particles;i++){
			for(j=0;j<n_cages;j++){
				Hopping_Rate_Now[i][j]=0;
				Hopping_Rate_Cumulative[i][j]=0;
			}
		}
		Hopping_Rate_Total=0;
		for(particle=0;particle<n_particles;particle++){
			cage=MoleculeDescription2[particle];
			for(j=0;j<n_cages;j++){
				if(CageConnectivities[cage][j]!=0){
					if(CageOccupancy[j]==0){
						r=CageTypes[cage];
						v=CageTypes[j];
						Hopping_Rate_Now[particle][j]=Hoprates[r][v];
						Hopping_Rate_Total+=Hopping_Rate_Now[particle][j];
					}
				}
			}
		}
		tau=-log(double((rand()%1048575)+1)/1048576.0)/Hopping_Rate_Total;
		a=double(((rand()%1048575)+1))/1048576.0;//a selects the event to execute
//		cout<<a<<endl;
		a*=Hopping_Rate_Total;
		Hopping_Rate_Cumulative[0][0]=Hopping_Rate_Now[0][0];//construct cumulative hopping rate matrix
		for(cage=1;cage<n_cages;cage++){
			Hopping_Rate_Cumulative[0][cage]=Hopping_Rate_Cumulative[0][cage-1]+Hopping_Rate_Now[0][cage];
		}
		if(n_particles>1){
			for(particle=1;particle<n_particles;particle++){
				Hopping_Rate_Cumulative[particle][0]=Hopping_Rate_Cumulative[particle-1][n_cages-1]+Hopping_Rate_Now[particle][0];
				for(cage=1;cage<n_cages;cage++){
					Hopping_Rate_Cumulative[particle][cage]=Hopping_Rate_Cumulative[particle][cage-1]+Hopping_Rate_Now[particle][cage];
				}
			}
		}
		done=0;
		upper_bound=n_cages*n_particles-1;
		lower_bound=0;
		while(done!=1){
			target=(upper_bound+lower_bound)/2;
			if (a<Hopping_Rate_Cumulative[target/n_cages][target%n_cages]){
				if(target==0)
			done=1;
				else if(a>Hopping_Rate_Cumulative[(target-1)/n_cages][(target-1)%n_cages])
					done=1;
				else
					upper_bound=target-1;
			}
			else
				lower_bound=target+1;
		}
		
		selected_particle=int(target/n_cages);
		selected_cage=target%n_cages;
		t=t+tau;
//		cout<<selected_particle<<" "<<selected_cage<<" "<<target<<endl;

		//----------update record----------
		if(t>0){
			while(t_record_max<t){
				record_index++;
				t_record_max+=t_increment;
				t_record[record_index]=t_record[record_index-1]+t_increment;
			//	cout << t_record_max/t_final*100<<" percent complete! "<<endl;
				for(particle=0;particle<n_particles;particle++){
					for(r=0;r<7;r++){
						xyz_record[particle][r][record_index]=MoleculeDescription[particle][r];
					}
				}
			}
		}
		//------------update reality---------
		CageOccupancy[MoleculeDescription2[selected_particle]]--;
		CageOccupancy[selected_cage]++;
		MoleculeDescription[selected_particle][0]=CageCenters[selected_cage][0];
		MoleculeDescription[selected_particle][1]=CageCenters[selected_cage][1];
		MoleculeDescription[selected_particle][2]=CageCenters[selected_cage][2];
		MoleculeDescription[selected_particle][3]+=double(CageConnectivities2[MoleculeDescription2[selected_particle]][selected_cage][0]);
		MoleculeDescription[selected_particle][4]+=double(CageConnectivities2[MoleculeDescription2[selected_particle]][selected_cage][1]);
		MoleculeDescription[selected_particle][5]+=double(CageConnectivities2[MoleculeDescription2[selected_particle]][selected_cage][2]);
		MoleculeDescription2[selected_particle]=selected_cage;
		
//following checks that Cage Occupancies are accounted correctly:
//		for(j=0;j<n_cages;j++){
//			done=CageOccupancy[j];
//			for(i=0;i<n_particles;i++){
//				if(MoleculeDescription2[i]==j){
//					done--;
//				}
//			}
//			if(done!=0){
//				cout<<"big problem"<<endl;
//			}
//		}
	}
	//----------KMC simulation over----------
	cout<<"KMC simulation complete"<<endl;
	//----------Now MSD analysis----------
	int t_index;
	double x_record[n_particles][int(t_final/t_increment)];
	double y_record[n_particles][int(t_final/t_increment)];
	double z_record[n_particles][int(t_final/t_increment)];
	for(j=0;j<int(t_final/t_increment);j++){
		for(particle=0;particle<n_particles;particle++){
			x_record[particle][j]=0;
			y_record[particle][j]=0;
			z_record[particle][j]=0;
		}
	}
	//unfold PBC:
//	cout<<int(t_final/t_increment)<<endl;
	for(particle=0;particle<n_particles;particle++){
//		cout<<"now recording particle: "<<particle<< " ";
		for(t_index=2;t_index<int(t_final/t_increment);t_index++){
			x_record[particle][t_index-1]=xyz_record[particle][0][t_index]+xyz_record[particle][3][t_index]*BoxDimensions[0];
			y_record[particle][t_index-1]=xyz_record[particle][1][t_index]+xyz_record[particle][4][t_index]*BoxDimensions[1];
			z_record[particle][t_index-1]=xyz_record[particle][2][t_index]+xyz_record[particle][5][t_index]*BoxDimensions[2];
//			cout<<x_record[particle][t_index-1]<<" ";
//			cout<<"y="<<y_record[particle][t_index-1]<<" ";
//			cout<<"z="<<z_record[particle][t_index-1]<<" ";
		}
//		cout<<"."<<endl;
	}
	//figure out which time lengths to investigate, logarithmically spaced by 2
	i=int(log2((t_final*0.95)/t_increment));
	double time_lengths_to_investigate[i];
	time_lengths_to_investigate[0]=t_increment;
	for(j=1;j<i;j++){
		time_lengths_to_investigate[j]=time_lengths_to_investigate[j-1]*2;
	}
	//first, self diffusion coefficient
	cout<<endl<<"MSD_self:"<<endl;
	int t_start_index,number_to_skip;
	double Square_Displacement_x, Square_Displacement_y, Square_Displacement_z;
	double MSD_self_x[i], MSD_self_y[i], MSD_self_z[i];
	for(t_index=0;t_index<i;t_index++){
		number_to_skip=int(time_lengths_to_investigate[t_index]/t_increment);
		z=int((t_final-time_lengths_to_investigate[t_index])/t_increment-2.0);
		if(z>250)
			z=250;
		Square_Displacement_x=Square_Displacement_y=Square_Displacement_z=0;
		j=0;
		for(t_start_index=1;t_start_index<z;t_start_index++){
			for(particle=0;particle<n_particles;particle++){
				a=(x_record[particle][t_start_index]-x_record[particle][t_start_index+number_to_skip]);
				Square_Displacement_x+=a*a;
				a=(y_record[particle][t_start_index]-y_record[particle][t_start_index+number_to_skip]);
				Square_Displacement_y+=a*a;
				a=(z_record[particle][t_start_index]-z_record[particle][t_start_index+number_to_skip]);
				Square_Displacement_z+=a*a;
				j++;
			}
		}
		MSD_self_x[t_index]=Square_Displacement_x/j;
		MSD_self_y[t_index]=Square_Displacement_y/j;
		MSD_self_z[t_index]=Square_Displacement_z/j;
		cout<< MSD_self_x[t_index]<< " "<<MSD_self_y[t_index]<< " "<<MSD_self_z[t_index]<<" "<<time_lengths_to_investigate[t_index]<<endl;
	}
	cout << endl<<"MSD_Collective:"<<endl;
	// Next, collective diffusion coefficient
	double x_record_com[int(t_final/t_increment)], y_record_com[int(t_final/t_increment)], z_record_com[int(t_final/t_increment)];
	for(t_index=0;t_index<int(t_final/t_increment);t_index++){
		x_record_com[t_index]=0.0;
		y_record_com[t_index]=0.0;
		z_record_com[t_index]=0.0;
		for(particle=0;particle<n_particles;particle++){
			x_record_com[t_index]+=x_record[particle][t_index];
			y_record_com[t_index]+=y_record[particle][t_index];
			z_record_com[t_index]+=z_record[particle][t_index];
		}
	}
	double MSD_collective_x[i], MSD_collective_y[i], MSD_collective_z[i];
	for(t_index=0;t_index<i;t_index++){
		number_to_skip=int(time_lengths_to_investigate[t_index]/t_increment);
		z=int((t_final-time_lengths_to_investigate[t_index])/t_increment-2.0);
		if(z>250)
			z=2500;
		j=0;
		Square_Displacement_x=Square_Displacement_y=Square_Displacement_z=0;
		for(t_start_index=1;t_start_index<z;t_start_index++){
			a=(x_record_com[t_start_index]-x_record_com[t_start_index+number_to_skip]);
			Square_Displacement_x+=a*a;
			a=(y_record_com[t_start_index]-y_record_com[t_start_index+number_to_skip]);
			Square_Displacement_y+=a*a;
			a=(z_record_com[t_start_index]-z_record_com[t_start_index+number_to_skip]);
			Square_Displacement_z+=a*a;
			j++;
		}
		MSD_collective_x[t_index]=Square_Displacement_x/j/n_particles;
		MSD_collective_y[t_index]=Square_Displacement_y/j/n_particles;
		MSD_collective_z[t_index]=Square_Displacement_z/j/n_particles;
		cout<< MSD_collective_x[t_index]<< " "<<MSD_collective_y[t_index]<< " "<<MSD_collective_z[t_index]<<" "<<time_lengths_to_investigate[t_index]<<endl;
	}
	double D_self[i+1];
	double D_collective[i+1];
	for(t_index=0;t_index<i;t_index++){
		D_self[t_index]=(MSD_self_x[t_index]+MSD_self_y[t_index]+MSD_self_z[t_index])/6/time_lengths_to_investigate[t_index];
		D_collective[t_index]=(MSD_collective_x[t_index]+MSD_collective_y[t_index]+MSD_collective_z[t_index])/6/time_lengths_to_investigate[t_index];
	}
	//First do self
	double D_self_final;
	double D_collective_final;
	done=0;
	double stdev;
	n=i;
	double sum_of_squares;
	sum_of_squares=0;
	double average_value;
	average_value=0;
	while(done==0){
		for(t_index=i-n;t_index<i;t_index++){	
			average_value+=D_self[t_index];
		}
		average_value=average_value/n;
		for(t_index=i-n;t_index<i;t_index++){	
			sum_of_squares+=(average_value-D_self[t_index])*(average_value-D_self[t_index]);
		}
		stdev=sqrt(sum_of_squares/n);
		if(stdev<average_value/10){
			done = 1;
		}
		else{
			if(n>=i-1){
				done = 1;
			}
			else{
				n++;
			}
		}
	}
	D_self_final=average_value;
	//then collective
	sum_of_squares=0;
	average_value=0;
	done = 0;
	while(done==0){
		for(t_index=i-n;t_index<i;t_index++){	
			average_value+=D_collective[t_index];
		}
		average_value=average_value/n;
//		cout <<"average value is "<<average_value<<endl;
		for(t_index=i-n;t_index<i;t_index++){	
			sum_of_squares+=(average_value-D_collective[t_index])*(average_value-D_collective[t_index]);
		}
		stdev=sqrt(sum_of_squares/n);
		if(stdev<average_value/10){
			done = 1;
		}
		else{
			if(n>=i-1){
				done = 1;
			}
			else{
				n++;
			}
		}

	}
	D_collective_final=average_value;
	
	D_self_global=D_self_final;
	D_collective_global=D_collective_final;
	return 0;
}
		
int main ()
{
	int i;
	int j;
	double D[5][2];
	for(i=0;i<5;i++){
		KMC();
		D[i][0]=D_self_global;
		D[i][1]=D_collective_global;
		cout << D[i][0]<<" "<<D[i][1]<<endl;
	}
	cout << "done with initial five studies"<<endl;
	double sum_of_values=0;
	for(i=0;i<5;i++){
		sum_of_values+=D[i][1];
		cout << D[i][0]<<" "<<D[i][1]<<endl;
	}
	double average_value;
	average_value=sum_of_values/5;
	double sum_of_squares=0;
	for(i=0;i<5;i++){
		sum_of_squares+=(average_value-D[i][1]);
	}
	double ste;
	ste=sqrt(sum_of_squares/5)/sqrt(5);
	double error_factor;
	error_factor=ste/average_value;
	double tolerance;
	tolerance=0.1;
	int number_of_iterations;
	if (error_factor<tolerance){
		number_of_iterations=5;
	}
	else{
		number_of_iterations=ceil((error_factor/tolerance)*(error_factor/tolerance));
	}
	double D2[number_of_iterations][2];
	for(i=0;i<number_of_iterations;i++){
		KMC();
		D2[i][0]=D_self_global;
		D2[i][1]=D_collective_global;
		cout << D2[i][0]<<" "<<D2[i][1]<<endl;
	}
	sum_of_values=0;
	for(i=0;i<number_of_iterations;i++){
		sum_of_values+=D2[i][0];
	}
	D_self_global=sum_of_values/number_of_iterations;
	sum_of_values=0;
	for(i=0;i<number_of_iterations;i++){
		sum_of_values+=D[i][1];
	}
	D_collective_global=sum_of_values/number_of_iterations;
	cout << endl<<D_self_global << " "<< D_collective_global<< endl;
	return 0;
}		
		
		
		
		
		
