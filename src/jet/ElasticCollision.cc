#include "ElasticCollision.h"
#include "scat.h"
#include <TRotation.h>
#include <TVector3.h>


std::uniform_real_distribution<> dis(0.0,1.0);
ElasticCollision::ElasticCollision(){};

void ElasticCollision::setter(double higher_temp0,double lower_temp0,double div_temp0,double max_energy0,double low_energy0,double energy_grid0){
    obj_low_temp=lower_temp0;
    obj_hig_temp=higher_temp0;
    obj_grid_temp=div_temp0;
    obj_hig_energy=max_energy0;
    obj_low_energy=low_energy0;
    obj_grid_energy=energy_grid0;
    iZ_Vector.SetXYZ(0,0,1);
    scattering_obj.setter(higher_temp0,lower_temp0,div_temp0,max_energy0,low_energy0,energy_grid0);
    scattering_obj.initialize_samplers(9999);
}
/*The MATTER module had pre-defination  of calculating the probability of scattring.*/

void ElasticCollision::elastic_kinematics(double temperature,int& pid0,int& pid2, int& pid3 ,double (&pc0)[4],double (&pc2)[4],double (&pc3)[4]){
	double r0,r1,r2,r3,r4,r5;
    double E1,E2,E3,E4,c23,c24,p1,p3,THETA2,THETA3,THETA4,PHI3,PHI4,PHI2,PHI23;
    double mc_sq;

    int parent_pid,hole_pid,daughter1_pid,daughter2_pid,pid_index;
    int energy_index, temp_index;
    int parton_type, proc; //heavy(1) or light(0)

    TVector3 P1;
    TVector3 P2;
    TVector3 P3;
    TVector3 P4;
    TRotation r;

    E1=pc0[0];
    energy_index=round((E1- obj_low_energy)/obj_grid_energy);
    temp_index=round((temperature-obj_low_temp)/obj_grid_temp);
    E1=scattering_obj.get_energy(energy_index);
    parent_pid=pid0;
    pid_index=floor(dis(generator)*6);
    if (parent_pid==21){
        parton_type=0;

		r0=scattering_obj.get_rate(temp_index,energy_index,6);
		r1=scattering_obj.get_rate(temp_index,energy_index,7);
		r2=scattering_obj.get_rate(temp_index,energy_index,8);  

		std::discrete_distribution<int> distribution0{r0,r1,r2};
		proc=distribution0(generator);
        if(proc==0){
            scattering_obj.get_sample(temp_index,energy_index,6,V);
            hole_pid=21;
            daughter1_pid=pid_list[pid_index];
            daughter2_pid=-pid_list[pid_index];
        } 
        else if(proc==1){
            scattering_obj.get_sample(temp_index,energy_index,7,V);
            hole_pid=21;
            daughter1_pid=21;
            daughter2_pid=-21;
        } 
		else{
            scattering_obj.get_sample(temp_index,energy_index,8,V);
            hole_pid=pid_list[pid_index];
            daughter1_pid=21;
            daughter2_pid=pid_list[pid_index];
        }            
    }
    else if(parent_pid>=1 and parent_pid<=3){
        parton_type=0;

		r0=scattering_obj.get_rate(temp_index,energy_index,0);
		r1=scattering_obj.get_rate(temp_index,energy_index,1);
		r2=scattering_obj.get_rate(temp_index,energy_index,2);
		r3=scattering_obj.get_rate(temp_index,energy_index,3);
		r4=scattering_obj.get_rate(temp_index,energy_index,4);
		r5=scattering_obj.get_rate(temp_index,energy_index,5);

		std::discrete_distribution<int> distribution0{r0,r1,r2,r3,r4,r5};
		proc=distribution0(generator);	    //not sure
        if(proc==0){
            scattering_obj.get_sample(temp_index,energy_index,0,V);
            hole_pid=-parent_pid;
            do {
                pid_index=floor(dis(generator)*6);
                daughter1_pid = pid_list[pid_index];
            } while (daughter1_pid == parent_pid);                
            daughter2_pid=-daughter1_pid;
        }
        else if(proc==1){
            scattering_obj.get_sample(temp_index,energy_index,1,V);
            hole_pid=-parent_pid;
            daughter1_pid=parent_pid;
            daughter2_pid=-parent_pid;                
        }
        else if(proc==2){
            scattering_obj.get_sample(temp_index,energy_index,2,V);
            hole_pid=parent_pid;
            daughter1_pid=parent_pid;
            daughter2_pid=parent_pid;
        } 
        else if(proc==3){
            scattering_obj.get_sample(temp_index,energy_index,3,V);
            hole_pid=-parent_pid;
            daughter1_pid=21;
            daughter2_pid=21;
        }
        else if(proc==4){
            scattering_obj.get_sample(temp_index,energy_index,4,V);
            hole_pid=21;
            daughter1_pid=parent_pid;
            daughter2_pid=21;
        }  
        else{ 
            scattering_obj.get_sample(temp_index,energy_index,5,V);
            daughter1_pid=parent_pid;
            do {
                pid_index=floor(dis(generator)*6);
                daughter2_pid = pid_list[pid_index];
            } while (daughter2_pid == parent_pid || daughter2_pid == -parent_pid);
            hole_pid=daughter2_pid;          
        }       
    }
     /*Heavy scattering part will be added here*/
    else if(parent_pid==4){
        parton_type=1;
		r0=scattering_obj.get_rate(temp_index,energy_index,9);
		r1=scattering_obj.get_rate(temp_index,energy_index,10);
		std::discrete_distribution<int> distribution0{r0,r1};
		proc=distribution0(generator);	
        if(proc==0){
            scattering_obj.get_sample(temp_index,energy_index,9,V);
            hole_pid=pid_list[pid_index];
            daughter1_pid=parent_pid;
            daughter2_pid=pid_list[pid_index];
        }
        else{
            scattering_obj.get_sample(temp_index,energy_index,10,V);
            hole_pid=21;
            daughter1_pid=parent_pid;
            daughter2_pid=21;
        }
        mc_sq=1.6129;
    }
    else{
        parton_type=1;
		r0=scattering_obj.get_rate(temp_index,energy_index,11);
		r1=scattering_obj.get_rate(temp_index,energy_index,12);
		std::discrete_distribution<int> distribution0{r0,r1};
		proc=distribution0(generator);	
        if(proc==0){
            scattering_obj.get_sample(temp_index,energy_index,11,V);
            scattering_obj.get_sample(temp_index,energy_index,9,V);
            hole_pid=pid_list[pid_index];
            daughter1_pid=parent_pid;
            daughter2_pid=pid_list[pid_index];
        }
        else{
            scattering_obj.get_sample(temp_index,energy_index,12,V);
            hole_pid=21;
            daughter1_pid=parent_pid;
            daughter2_pid=21;
        }
        mc_sq=17.4724;
    }

    if (parton_type==0){
        THETA2=V[0];
        THETA3=V[1];
        PHI23=V[2];
        E3=V[3];

		c23=cos(THETA2)*cos(THETA3)+sin(THETA2)*sin(THETA3)*cos(PHI23);
		E2=(E1*E3*(1-cos(THETA3)))/(E1*(1-cos(THETA2))-E3*(1-c23));
		E4=E1+E2-E3;
		THETA3=THETA3;
		PHI2=dis(generator);
		PHI3=PHI2-PHI23;
		THETA4= acos((E1+E2*cos(THETA2)-E3*cos(THETA3))/E4);
		PHI4=acos((E2*sin(THETA2)*cos(PHI2)-E3*sin(THETA3)*cos(PHI3))/(E4*sin(THETA4)));

        P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        P2.SetXYZ(E2*sin(THETA2)*cos(PHI2),E2*sin(THETA2)*sin(PHI2),E2*cos(THETA2));
		P3.SetXYZ(E3*sin(THETA3)*cos(PHI3),E3*sin(THETA3)*sin(PHI3),E3*cos(THETA3));
		P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));

		double s0=P1.Angle(iZ_Vector);
		TVector3 s1;
		s1=P1.Cross(iZ_Vector);

		r.Rotate(s0,s1);
        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
    }
    else{
		E2=V[0];
        THETA2=V[1];
        THETA4=V[2];
        PHI4=V[3];

        c24=sin(THETA2)*sin(THETA4)*cos(PHI4)+cos(THETA2)*cos(THETA4);
        p1=sqrt(E1*E1-mc_sq);
        E4=(E1*E2-p1*E2*cos(THETA2))/(E1-p1*cos(THETA4)+E2-E2*c24);
        E3=E1+E2-E4;
        p3=sqrt(E3*E3-mc_sq);
        THETA3=acos((p1+E2*cos(THETA2)-E4*cos(THETA4))/p3);
        PHI3=asin(-(E4*sin(PHI4)*sin(THETA4))/(p3*sin(THETA3)));
        PHI2=dis(generator);     
        r.RotateZ(PHI2);

        P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        P2.SetXYZ(E2*sin(THETA2),0,E2*cos(THETA2));
        P3.SetXYZ(p3*sin(THETA3)*cos(PHI3),p3*sin(THETA3)*sin(PHI3),p3*cos(THETA3));
        P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));

        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
    }

    double s0=P1.Angle(iZ_Vector);
    TVector3 s1;
    s1=P1.Cross(iZ_Vector);
    TRotation w0;
    w0.Rotate(s0,s1);
    
    P2=w0*P2;
    P3=w0*P3;
    P4=w0*P4;

    pc3[1]=E2;
    pc3[1]=P2.x();
    pc3[2]=P2.y();
    pc3[3]=P2.z();
    pid3=hole_pid;

    pc0[1]=E3;
    pc0[1]=P3.x();
    pc0[2]=P3.y();
    pc0[3]=P3.z();
    pid0=daughter1_pid;

    pc2[1]=E4;
    pc2[1]=P4.x();
    pc2[2]=P4.y();
    pc2[3]=P4.z();   
    pid2=daughter2_pid; 
}


