subroutine potential(eng,eps_cl_hcl,sig_cl_hcl,eps_cl_v,sig_cl_v,eps_hyd_hydron_cl,sig_hyd_hydron_cl,aof,bof,&
& a1,b1,cons1,cons2,repa,repb,repd,repc,eps_hhcl_oxw,sig_hhcl_oxw,eps_hhcl_hydw,sig_hhcl_hydw,eps_clhcl_hydw,&
&sig_clhcl_hydw,q_hhcl,b_clo_1,b_clo_2,b_clo_3,doo_clo,fac,alpha,fac2,beta,fac3,alph2,r1sw,r2sw) 
use hcl_mod
implicit none

!!! ===== Variables =====
real*8, allocatable, dimension(:) ::d,eig_vec  !(eigenvalue)
real*8, allocatable, dimension(:,:) ::v  !(eigen vector)
integer :: nrot,grstate,ground_state
real*8 :: r1sw,r2sw
real*8 :: term_hydron_cl_rep,q1clo,q2clo,q3clo,b_clo_1,b_clo_2,b_clo_3,doo_clo
real*8 :: cons1,cons2,temp_2
real*8 :: term_1, term_2 ,semper, min,eng,aof,bof,alpha,beta,fac,fac2,fac3,alph2,a1,b1, repa,repb,repd,repc 
real*8 :: sig_hhcl_oxw,eps_hhcl_oxw, sig_hhcl_hydw,eps_hhcl_hydw, eps_clhcl_hydw,sig_clhcl_hydw
real*8 :: sig_hhcl_oxw_12,sig_hhcl_oxw_6,eps_hhcl_oxw_4,sig_hhcl_hydw_12, sig_hhcl_hydw_6, eps_hhcl_hydw_4
real*8 :: eps_clhcl_hydw_4,sig_clhcl_hydw_12,sig_clhcl_hydw_6 
real*8 :: hydron_roh, hydron_ang, eps_oxh_oxw_4, sig_oxh_oxw_12, sig_oxh_oxw_6, recip_r_oxh_oxw
real*8 :: recip_r_hydrh_oxw, eps_hydrh_oxw_4, sig_hydrh_oxw_12, sig_hydrh_oxw_6, qhyd_o, qhyd_h
real*8 :: qwat_o, qwat_h, recip_r_oxh_wath, recip_r_hydh_wath, r_oxh_oxw 
real*8 :: hydron_hyd_wat_ox_dis, roh1_bond, roh2_bond,rhcl_0,hcl_3_body,hcl_oxw1,hcl_oxw2
real*8 :: temp, term_3, roh1_ang_roh2, eps_oxw_oxw_4, sig_oxw_oxw_12, sig_oxw_oxw_6, recip_r_oxw_oxw
real*8 :: recip_r_hydw_oxw, recip_r_hydw_hydw, sig_cl_w, eps_cl_w, sig_ox_wat, eps_ox_wat, eps_cl_oxw
real*8 :: sig_cl_oxw, r_cl_oxw, recip_r_cl_oxw, term_4, q_cl, recip_r_cl_hydw, sig_ox_hydron
real*8 :: eps_ox_hydron, term_5, sig_ox_hydron_cl, eps_ox_hydron_cl, r_cl_ox_hydron, recip_r_cl_ox_hydron       
real*8 :: qhyd_o_w, qhyd_h_w, r_cl_hyd_hydron, hydron_cl_pot, term_6, sig_cl_hcl, eps_cl_hcl         
real*8 :: oxw_oxw_recip, hydw_oxw_recip, hydw_hydw_recip, clhcl_oxw,clhcl_oxw_recip  
real*8 :: term_7, term_8, q_clhcl, clhcl_hydw_recip, q_hhcl, hhcl_oxw_recip, hhcl_hydw_recip           
real*8 :: term_9, h_cl_r, wat_hcl_pot, eps_hydh, sig_hydh, sig_hyd_hydron_cl, sig_hyd_hydron_cl_12
real*8 :: sig_hyd_hydron_cl_6, eps_hyd_hydron_cl, eps_hyd_hydron_cl_4    , term_51,term_52 
real*8 :: q_cl_v, sig_cl_v, eps_cl_v, eps_oxw, sig_oxw, eps_cl_oxw_4, sig_cl_oxw_6, sig_cl_oxw_12  
real*8 :: eps_oxh, sig_oxh,sig_ox_hydron_cl_12, sig_ox_hydron_cl_6, eps_ox_hydron_cl_4   
real*8 :: lj,coulomb,rep_o,rep_ho, qx,qy,qz,q1_sq,q2_sq,q3_sq, term_3_body   
real*8 :: sig_hyd_wat_cl,sig_hyd_wat_cl_12,sig_hyd_wat_cl_6,eps_hyd_wat_cl,eps_hyd_wat_cl_4
real*8 :: a1_sw,a2_sw,angoxclox,angohydrcl, swang, swr1, swr2,hhcl_oxw,term_vrep 
real*8 :: eps_oxh_hydw, eps_oxh_hydw_4, sig_oxh_hydw, sig_oxh_hydw_6, sig_oxh_hydw_12
integer :: i,j,k,l,m
integer :: hydron_ox, hydron_h, h_a, h_b, wat_h, roh1_index, roh2_index, roh3_index, wat_h2
integer :: wat_h3, test_hyd, test_hyd_2, test_hyd_3, h_hcl, wat_h4   

!!! ===== Variables Off-Diagonal =====
real*8 :: term_10,term_11,term_12, recip_zund_ox_1_oxw, recip_zund_ox_2_oxw, qzund_o, qzund_h,qzund_h_s
real*8 :: r_sharz_hyd_oxw, r_z_ox_h_wat_1, r_z_ox_h_wat_2, r_zsh_wath, rzh_oxw_1, rzh_oxw_2, rzh_hydw
real*8 :: roo_dist, scale_fact, v_ij_ex, off_hyd_hyd_pot, v_ij_const,hcl_h_zund_dist,off_hyd_hcl_pot
real*8 :: q, q1,q2,q3
integer :: counter, zund_ox_1,zund_ox_2,zund_h_shar, hydron_hyd_1, hydron_hyd_2,zund_hcl_hyd

!!! ===== Temp =====
real*8 :: eps_cl_oxw_w,eps_cl_oxw_4_w,sig_cl_oxw_w, sig_cl_oxw_12_w,sig_cl_oxw_6_w

!!! ===== MS-EVB 3.2 Parameters =====
eps_oxh=0.098609686d0	! kcal/mol
sig_oxh=3.118508d0	! angstroms
eps_hydh=0.000040458d0	!kcal/mol
sig_hydh=0.0d0		! angstroms
eps_oxh_hydw=3.d0	!kcal/mol
sig_oxh_hydw=1.6d0	!angstroms
eps_oxh_hydw_4=eps_oxh_hydw*4.d0
sig_oxh_hydw_6=sig_oxh_hydw**6.d0
sig_oxh_hydw_12=sig_oxh_hydw**12.d0
qhyd_o=-0.32d0
qhyd_h=0.44d0
!!! ===== off diagonal parameters =====
!!!!! ===== MSEVB3.2 PARAMETERS =====
v_ij_const=-21.064268d0		!kcal/mol
qzund_o=-0.0895456d0		! e 
qzund_h=0.0252683d0		!e
qzund_h_s=0.0780180d0		! e

!!! ===== Voth Cl Parameters =====
q_cl_v=-1.0d0 ! q, partial charge 

!!! ===== SPC/FW Parameters =====
qwat_o=-0.82d0
qwat_h=0.41d0
eps_oxw=0.1554235d0	!kcal/mol 
sig_oxw=3.165492d0 !Angstroms
eps_oxw_oxw_4=(0.15542353d0)*4.d0	!kcal/mol
sig_oxw_oxw_12=(3.165492d0)**12d0		!angstroms to the 12
sig_oxw_oxw_6=(3.165492d0)**6d0		!angstroms to the 6

!!! ===== Wick Parameters =====
!q_clhcl=-0.23d0		! e 
!q_hhcl=0.23d0		! e
q_clhcl=-q_hhcl

!!! ===== Lennard-Jones Paramters/Combining Rules/Raised to the Correct Power =====
!!! new may 1st h of hcl to oxygen of water 
eps_hhcl_oxw_4=eps_hhcl_oxw*4.d0
sig_hhcl_oxw_12=sig_hhcl_oxw**12.d0
sig_hhcl_oxw_6=sig_hhcl_oxw**6.d0
!!! new may 1st h of hcl to hydrogen of water
eps_hhcl_hydw_4=eps_hhcl_hydw*4.d0
sig_hhcl_hydw_12=sig_hhcl_hydw**12.d0
sig_hhcl_hydw_6=sig_hhcl_hydw**6.d0
!!! new may 4th cl of hcl to hydrogen of water
eps_clhcl_hydw_4=eps_clhcl_hydw*4.d0
sig_clhcl_hydw_12=sig_clhcl_hydw**12.d0
sig_clhcl_hydw_6=sig_clhcl_hydw**6.d0

 
!!! ##### oxygen of hydronium to oxygen of water ####
eps_oxh_oxw_4=(sqrt(eps_oxh*eps_oxw))*4.d0
sig_oxh_oxw_12=((sig_oxh+sig_oxw)/2.d0)**12.d0 
sig_oxh_oxw_6=((sig_oxh+sig_oxw)/2.d0)**6.d0 

!!! ##### hydrogen of hydronium to oxygen of water #####
eps_hydrh_oxw_4=(sqrt(eps_hydh*eps_oxw))*4.d0
sig_hydrh_oxw_12=((sig_hydh+sig_oxw)/2.d0)**12.d0
sig_hydrh_oxw_6=((sig_hydh+sig_oxw)/2.d0)**6.d0 

!!! ##### cl- and water oxygen mixing rules #####
eps_cl_oxw=sqrt(eps_cl_v*eps_oxw)
eps_cl_oxw_4=eps_cl_oxw*4.d0
sig_cl_oxw=((sig_cl_v+sig_oxw)/2.d0)
sig_cl_oxw_12=sig_cl_oxw**12.d0 
sig_cl_oxw_6=sig_cl_oxw**6.d0 

!!! ##### Cl of HCl (wick) to oxygen of water #####
eps_cl_oxw_w=sqrt(eps_cl_hcl*eps_oxw)
eps_cl_oxw_4_w=eps_cl_oxw_w*4.d0
sig_cl_oxw_w=((sig_cl_hcl+sig_oxw)/2.d0)
sig_cl_oxw_12_w=sig_cl_oxw_w**12.d0 
sig_cl_oxw_6_w=sig_cl_oxw_w**6.d0 

!!! ##### cl- and oxygen of hydronium ####
sig_ox_hydron_cl=((sig_oxh+sig_cl_v)/2.d0)
eps_ox_hydron_cl=sqrt(eps_oxh*eps_cl_v)
sig_ox_hydron_cl_12=sig_ox_hydron_cl**12.d0 
sig_ox_hydron_cl_6=sig_ox_hydron_cl**6.d0 
eps_ox_hydron_cl_4=eps_ox_hydron_cl*4.d0

!!! ##### cl- and hydrogen of hydronium #####
sig_hyd_hydron_cl_12=sig_hyd_hydron_cl**12.d0 
sig_hyd_hydron_cl_6=sig_hyd_hydron_cl**6.d0 
eps_hyd_hydron_cl_4=eps_hyd_hydron_cl*4.d0

!!! ===== LJ Parameters cl- to h of water =====
!!! ===== J. Chem. Phys. 122, 144105, 2005 =====
!!! ===== DOI:10.1063/1.1881092 =====
sig_hyd_wat_cl=1.d0 ! angstroms
sig_hyd_wat_cl_12=sig_hyd_wat_cl**12.d0 
sig_hyd_wat_cl_6=sig_hyd_wat_cl**6.d0 
eps_hyd_wat_cl=0.0001d0 !kcal/mol
eps_hyd_wat_cl_4=eps_hyd_wat_cl*4.d0


!!! fitted two body parameters 
!eps_cl_hcl=0.99979D+00 
!sig_cl_hcl=0.40621D+01 
!eps_cl_v= 0.78156D-01 
!sig_cl_v= 0.16410D+01 
!eps_hyd_hydron_cl= 0.18677D-01 
!sig_hyd_hydron_cl=0.31263D+00 
!aof=-0.20309716D+02 
!bof=0.93201D+00 
!a1=-0.15855D+02 
!b1=0.52438D+01  
!cons1=0.53331D+01 
!cons2=0.10949D+02 
!repa=0.31314D+02 
!repb=0.41656653D-02 
!repd=0.13777D+01 
!repc=0.18784D+01 
!eps_hhcl_oxw=0.87460D-01 
!sig_hhcl_oxw= 0.18961D+00  
!eps_hhcl_hydw= 0.98917D+00  
!sig_hhcl_hydw= 0.16551D+01 
!eps_clhcl_hydw=0.11328550D+00 
!sig_clhcl_hydw=0.28580D+00 
!q_hhcl=0.43907D+00 
!b_clo_1=0.17476D+02 
!b_clo_2=0.10398D+00 
!b_clo_3=0.23576D+01 
!doo_clo=0.22587D+01 
!q_clhcl=-q_hhcl

!!!HCl Hydronium States
! term_1 is the intramolecular energy for the hydronium
! term_2 is the intermolecular interactions between water and hydronium
! term_3 is spc/fw for water-water intra and intermoleculer interactions for hydronium cl states
! term_4 is cl ion vander waals interaction and coulomb interactions with water
! term_5 is chlorine ion and hydronium vander waals and coulomb interaction
! term_3_body is 3 body cl- to oxygen
!!! HCL water states
! term_6 is water spc/fw for water hcl states
! term_7 is hcl water lj 
! term_8 is coulomb interactions->cl of hcl to water and h of hcl to water
! term_9 is the hcl morse potential 
! term_vrep is hcl to oxygen repulsion
do i=1, num_states
!print*, "state ",i
	if (i .LE. hydronium_states) then ! hydronium, Cl-, and water
		term_1=0.d0
		term_2=0.d0
		term_3=0.d0	
		term_4=0.d0
		term_5=0.d0
		term_51=0.d0
		term_52=0.d0
		term_3_body=0.d0
		temp=0.d0
		lj=0.d0
		coulomb=0.d0
		rep_o=0.d0
		rep_ho=0.d0
		term_hydron_cl_rep=0.d0
		hydron_ox=hydronium_oxygen(i) !sets the hydronium oxygen for those states 
		!!! ===== Hydronium Intramolecular =====
		do j=1,3
			hydron_roh=ox_hyd_bond(hydron_ox,hydron_con(hydron_ox,3*i-3+j))
			term_1=term_1+v_hydron_intra_bond(hydron_roh)
			!!! checked may 10 2017	
			!!! Calculates all of the hydronium bond intramolecular reactions
			do k=j+1,3
				h_a=hydron_con(hydron_ox, 3*i-3+j) ! find the index of the first hydrogen 
				h_b=hydron_con(hydron_ox, 3*i-3+k) ! finds the index of the second hydrogen next over 
				hydron_ang=hyd_ox_hyd_ang(hydron_ox,h_a,h_b) !find the angle between the two 
				term_1=term_1+v_hydron_intra_ang(hydron_ang)
				!!! checked may 10 2017	
			end do
		end do 
		!print*, "term_1", term_1		
		!!! ===== Hydronium-Water Intermolecular =====
		do j=1,num_ox ! sums over the water molecules
			if (j .NE. hydron_ox) then ! this is the amount of waters when hydronium present, takes out hydronium oxygen
				r_oxh_oxw=ox_ox_dis(hydron_ox,j)
				recip_r_oxh_oxw=recip_ox_ox_dis(hydron_ox, j) !!! distance between oxygen of hydronium and oxygen of water
				term_2=term_2+v_inter_lj(recip_r_oxh_oxw,sig_oxh_oxw_12, sig_oxh_oxw_6, eps_oxh_oxw_4) ! LJ for oxh to ox wat
				!!! Checked may 10 2017 	
				term_2=term_2+v_inter_coulomb(recip_r_oxh_oxw, qhyd_o, qwat_o)
				!!! Checked may 10 2017 
				do k=1,3 !!! Sums over the hydrogens connected to the water oxygen
					wat_h=hydron_con(j,3*i-3+k)
					if (wat_h .NE. 0) then
						recip_r_oxh_wath=recip_ox_hyd_bond(hydron_ox, wat_h)
						term_2=term_2+v_inter_coulomb(recip_r_oxh_wath, qhyd_o, qwat_h) 
						!!! checked may 10 2017	
						term_2=term_2+v_inter_lj(recip_r_oxh_wath, sig_oxh_hydw_12, sig_oxh_hydw_6, eps_oxh_hydw_4)
						!!! checked may 10 2017	
					end if
				end do
				do k=1,3 !!! Sums over the hydrogens of hydronium
					hydron_h=hydron_con(hydron_ox,3*i-3+k)
					recip_r_hydrh_oxw=recip_ox_hyd_bond(j,hydron_h) !finds recip bond dis from hydronium h to water oxygens
					term_2=term_2+v_inter_lj(recip_r_hydrh_oxw, sig_hydrh_oxw_12, sig_hydrh_oxw_6, eps_hydrh_oxw_4) ! LJ for hydh to ox wat
					! checked may 10 2017	
					term_2=term_2+v_inter_coulomb(recip_r_hydrh_oxw,qhyd_h, qwat_o)
					!!! Checked may 10 2017					
					hydron_hyd_wat_ox_dis=ox_hyd_bond(j,hydron_h)
					term_2=term_2+v_inter_hydo_rep_msevb3(hydron_hyd_wat_ox_dis) !ADD SWITCHING FUNCTION 
					!!! Checked may 10 2017	
					if (k == 1 ) then
						qx=((xx_O(hydron_ox)+xx_O(j))/2.d0)-xx_H(hydron_h)
						qy=((yy_O(hydron_ox)+yy_O(j))/2.d0)-yy_H(hydron_h)
						qz=((zz_O(hydron_ox)+zz_O(j))/2.d0)-zz_H(hydron_h)
						q1_sq=(qx**2)+(qy**2)+(qz**2) 	
					else if (k == 2) then
						qx=((xx_O(hydron_ox)+xx_O(j))/2.d0)-xx_H(hydron_h)
						qy=((yy_O(hydron_ox)+yy_O(j))/2.d0)-yy_H(hydron_h)
						qz=((zz_O(hydron_ox)+zz_O(j))/2.d0)-zz_H(hydron_h)
						q2_sq=(qx**2)+(qy**2)+(qz**2)   
					else if (k == 3) then
						qx=((xx_O(hydron_ox)+xx_O(j))/2.d0)-xx_H(hydron_h)
						qy=((yy_O(hydron_ox)+yy_O(j))/2.d0)-yy_H(hydron_h)
						qz=((zz_O(hydron_ox)+zz_O(j))/2.d0)-zz_H(hydron_h)
						q3_sq=(qx**2)+(qy**2)+(qz**2)	
					end if
					do l=1,3
						wat_h=hydron_con(j,3*i-3+l)
						if (wat_h .NE. 0) then !!! 
							recip_r_hydh_wath=recip_hyd_hyd_dist(hydron_h,wat_h)
							term_2=term_2+v_inter_coulomb(recip_r_hydh_wath, qhyd_h, qwat_h) 
      						!!! Checked may 11 2017							
						end if
					end do
				end do
				term_2=term_2+v_inter_oo_rep_msevb3(r_oxh_oxw,q1_sq,q2_sq,q3_sq) !ADD SWITCHING FUNCTION
				!!! Checked may 10 may 2017		
				!!! ===== Water only potential ===== !!!
				!!! ===== Intramolecular ===== !!!
				roh1_index=hydron_con(j,3*i-2) !finds all possible water bonds
				roh2_index=hydron_con(j,3*i-1)
				roh3_index=hydron_con(j,3*i)		
				!!! if statements essentially work around the fact that one index should be zero
				if (roh1_index == 0) then
					roh1_bond=ox_hyd_bond(j,roh2_index)
					roh2_bond=ox_hyd_bond(j,roh3_index)
					roh1_ang_roh2=hyd_ox_hyd_ang(j,roh2_index,roh3_index)
					term_3=term_3+v_intra_spcfw(roh1_bond,roh2_bond,roh1_ang_roh2)
				else if (roh2_index == 0) then
					roh1_bond=ox_hyd_bond(j,roh1_index)
					roh2_bond=ox_hyd_bond(j,roh3_index) 
					roh1_ang_roh2=hyd_ox_hyd_ang(j,roh1_index,roh3_index)
					term_3=term_3+v_intra_spcfw(roh1_bond,roh2_bond,roh1_ang_roh2)
				else if (roh3_index == 0) then
					roh1_bond=ox_hyd_bond(j,roh1_index)
					roh2_bond=ox_hyd_bond(j,roh2_index)
					roh1_ang_roh2=hyd_ox_hyd_ang(j, roh1_index, roh2_index)
					term_3=term_3+v_intra_spcfw(roh1_bond,roh2_bond,roh1_ang_roh2)
				end if	
				!!! Checked may 11 2017	
				!!! ===== Intermolecular --> LJ and Coulomb Potential ===== !!!	
				do k=j+1,num_ox
					if (k .NE. hydron_ox) then
						recip_r_oxw_oxw=recip_ox_ox_dis(j,k)
						term_3=term_3+v_inter_lj(recip_r_oxw_oxw, sig_oxw_oxw_12, sig_oxw_oxw_6, eps_oxw_oxw_4)
						!!! checked may 11 2017	
						term_3=term_3+v_inter_coulomb(recip_r_oxw_oxw,qwat_o,qwat_o)
						!!! Checked may 11 2017	
						do l=1,3
							wat_h2=hydron_con(j,3*i-3+l)
							if (wat_h2 .NE. 0) then
								recip_r_hydw_oxw=recip_ox_hyd_bond(k,wat_h2)
								term_3=term_3+v_inter_coulomb(recip_r_hydw_oxw,qwat_o,qwat_h)
								!!! Checked may 12 2017	
							end if
						end do							
						do l=1,3
							wat_h3=hydron_con(k,3*i-3+l)
							if (wat_h3 .NE. 0) then
								recip_r_hydw_oxw=recip_ox_hyd_bond(j,wat_h3)
								term_3=term_3+v_inter_coulomb(recip_r_hydw_oxw,qwat_o,qwat_h)
								!!! Checked may 12 2017	
								do m=1,3
									wat_h4=hydron_con(j,3*i-3+m)
									if (wat_h4 .NE. 0) then
										recip_r_hydw_hydw=recip_hyd_hyd_dist(wat_h3,wat_h4)
										term_3=term_3+v_inter_coulomb(recip_r_hydw_hydw,qwat_h,qwat_h)
										!!! checked may 12 2017	
									end if
								end do
							end if	
						end do	
						!!! 3_body ox_wat-> cl -> oxwat
						swr1=cl_o_dis(j)
						swr2=cl_o_dis(k) 
						term_3_body=term_3_body+rk3(swr1,swr2,fac,alpha,alpha)
						!!!Checked July 10 2017	
					end if
				end do
				!!! ===== chlorine anion water interaction =====
				r_cl_oxw=cl_o_dis(j)
				recip_r_cl_oxw=recip_cl_o_dis(j) 	
				term_4=term_4+v_inter_lj(recip_r_cl_oxw,sig_cl_oxw_12,sig_cl_oxw_6,eps_cl_oxw_4)			
				!!! Checked may 13 2017	
				term_4=term_4+v_inter_coulomb(recip_r_cl_oxw,qwat_o,q_cl_v)
				!!!Checked may 13 2017	
				do k=1,3
					if (hydron_con(j,3*i-3+k) .NE. 0) then
						recip_r_cl_hydw=recip_cl_h_dis(hydron_con(j,3*i-3+k))
						term_4=term_4+v_inter_coulomb(recip_r_cl_hydw,qwat_h,q_cl_v)
						!!! Checked may 13 2017	
						term_4=term_4+v_inter_lj(recip_r_cl_hydw,sig_hyd_wat_cl_12,sig_hyd_wat_cl_6,eps_hyd_wat_cl_4)	
						!!! Checked may 13 2017	
					end if
				end do

			end if	! ends if going through the oxygens that are not the hydronium oxygen
			end do	! end j through the total number of water molecules  
			!!! ===== chlorine anion hydronium interaction =====
			r_cl_ox_hydron=cl_o_dis(hydron_ox)
			recip_r_cl_ox_hydron=recip_cl_o_dis(hydron_ox)
			term_5=term_5+v_inter_lj(recip_r_cl_ox_hydron, sig_ox_hydron_cl_12, sig_ox_hydron_cl_6, eps_ox_hydron_cl_4)	
			!!! Checked May 14 2017	
			term_5=term_5+v_inter_coulomb(recip_r_cl_ox_hydron,qhyd_o,q_cl_v)
			!!! Checked May 14 2017	
			!!! below rep added july 3rd
                        do k=1,3 !!! Sums over the hydrogens of hydronium
				hydron_h=hydron_con(hydron_ox,3*i-3+k)
                                        if (k == 1 ) then
                                                qx=((xx_O(hydron_ox)+xx_cl)/2.d0)-xx_H(hydron_h)
                                                qy=((yy_O(hydron_ox)+yy_cl)/2.d0)-yy_H(hydron_h)
                                                qz=((zz_O(hydron_ox)+zz_cl)/2.d0)-zz_H(hydron_h)
                                                q1clo=(qx**2)+(qy**2)+(qz**2)
                                        else if (k == 2) then
                                                qx=((xx_O(hydron_ox)+xx_cl)/2.d0)-xx_H(hydron_h)
                                                qy=((yy_O(hydron_ox)+yy_cl)/2.d0)-yy_H(hydron_h)
                                                qz=((zz_O(hydron_ox)+zz_cl)/2.d0)-zz_H(hydron_h)
                                                q2clo=(qx**2)+(qy**2)+(qz**2)
                                        else if (k == 3) then
                                                qx=((xx_O(hydron_ox)+xx_cl)/2.d0)-xx_H(hydron_h)
                                                qy=((yy_O(hydron_ox)+yy_cl)/2.d0)-yy_H(hydron_h)
                                                qz=((zz_O(hydron_ox)+zz_cl)/2.d0)-zz_H(hydron_h)
                                                q3clo=(qx**2)+(qy**2)+(qz**2)
                     	               end if
			
			end do	
			term_hydron_cl_rep=v_inter_rep_hydron_cl(r_cl_ox_hydron,q1clo,q2clo,q3clo,b_clo_1,b_clo_2,b_clo_3,doo_clo)	
			!!! Checked July 10th 2017 
			do j=1,3
				r_cl_hyd_hydron=recip_cl_h_dis(hydron_con(hydron_ox,3*i-3+j))
				term_5=term_5+v_inter_lj(r_cl_hyd_hydron, sig_hyd_hydron_cl_12, sig_hyd_hydron_cl_6,eps_hyd_hydron_cl_4)	
				!!! Checked May 14 2017	
				term_5=term_5+v_inter_coulomb(r_cl_hyd_hydron,qhyd_h,q_cl_v)
				!!! Checked May 14 2017	
			end do	
			do j=1,num_ox
				if (j .ne. hydron_ox) then
				angohydrcl=1.d0	
				swr1=cl_o_dis(j)
				swr2=cl_o_dis(hydron_ox)
				swang=ox_cl_ox(hydron_ox,j)
				term_3_body=term_3_body+rk3_hyd(swr1,swr2,fac2,alpha,beta)
				!!! Checked July 10 2017	
				end if
			end do       
  !                 hydron_cl_pot=0.d0
			hydron_cl_pot=term_1+term_2+term_3+term_4+term_5 + term_3_body+term_hydron_cl_rep 
		hamil_matrix(i,i)=hydron_cl_pot 
		!print*, "term hydron cl rep", term_hydron_cl_rep
		!print*, "term_1", term_1 !this is good 
		!print*, "term_2", term_2 !this is good 
		!print*, "term_3", term_3  !this is good 
		!print*, "term_3_body",term_3_body  !This is good 
		!print*, "term_4", term_4 ! this is good 
		!print*, "term_5", term_5	! this is good 
		!print*, "total potential hydronium state", hydron_cl_pot ! this is good    
	else if (i .GT. hydronium_states) then ! HCl, water
	term_6=0.d0
	term_7=0.d0
	term_8=0.d0
	hcl_3_body=0.d0
	term_vrep=0.d0
	temp=0.d0	
	!!! ===== HCl intramolecular ===== !!! 
		h_hcl=wat_con(num_ox+1,(3*i)-(3*hydronium_states)-2)
		if (h_hcl == 0) then
			print*, "h of hcl is not defined"
                        stop
                end if
		h_cl_r=cl_h_dis(h_hcl)
		term_9=v_intra_hcl_morse(h_cl_r) 
		!!! checked june 7th 
		!!! Water only potential intramolecular ===== !!!
		do j=1,num_ox
			roh1_index=wat_con(j,(3*i)-(3*hydronium_states)-3+1)
			roh2_index=wat_con(j,(3*i)-(3*hydronium_states)-3+2)
			roh3_index=wat_con(j,(3*i)-(3*hydronium_states)-3+3)
			if (roh1_index == 0) then 
				roh1_bond=ox_hyd_bond(j,roh2_index)
				roh2_bond=ox_hyd_bond(j,roh3_index)
				roh1_ang_roh2 = hyd_ox_hyd_ang(j,roh2_index,roh3_index) 
			else if (roh2_index == 0) then
				roh1_bond=ox_hyd_bond(j,roh1_index)
				roh2_bond=ox_hyd_bond(j,roh3_index) 
				roh1_ang_roh2=hyd_ox_hyd_ang(j,roh1_index,roh3_index)  
			else if (roh3_index == 0) then
				roh1_bond=ox_hyd_bond(j,roh1_index)
				roh2_bond=ox_hyd_bond(j,roh2_index)
				roh1_ang_roh2=hyd_ox_hyd_ang(j,roh1_index,roh2_index)    
			end if 
			term_6=term_6+v_intra_spcfw(roh1_bond,roh2_bond,roh1_ang_roh2)
			!!! Checked may 14 2017	
			!!! ===== Intermolecular water terms ===== !!!	
			do k=j+1, num_ox
				oxw_oxw_recip=recip_ox_ox_dis(j,k)
				!!! Following gets all of the oxygen of j to oxygen of k LJ and Coulomb
				term_6=term_6+v_inter_lj(oxw_oxw_recip,sig_oxw_oxw_12,sig_oxw_oxw_6,eps_oxw_oxw_4)
				!!! Checked may 14 2017	
				term_6=term_6+v_inter_coulomb(oxw_oxw_recip,qwat_o,qwat_o)
				!!! Checked may 14 2017	
				do l=1,3
					test_hyd=wat_con(j,(3*i)-(3*hydronium_states)-3+l)
					if (test_hyd .NE. 0) then
						hydw_oxw_recip=recip_ox_hyd_bond(k,test_hyd)
						term_6=term_6+v_inter_coulomb(hydw_oxw_recip,qwat_o,qwat_h)
						!!!Checked may 14 2017		
				end if
				end do	
				do l=1,3
					test_hyd=wat_con(k,(3*i)-(3*hydronium_states)-3+l)
					if (test_hyd .NE. 0) then
						hydw_oxw_recip=recip_ox_hyd_bond(j,test_hyd)
						!!! The following should get all of the oxygen of j to the hydrogens
						!!! From oxygen k
						term_6=term_6+v_inter_coulomb(hydw_oxw_recip,qwat_o,qwat_h)
						!!!Checked may 14 2017		
						do m=1,3
							test_hyd_2=wat_con(j,(3*i)-(3*hydronium_states)-3+m)
							if (test_hyd_2 .NE. 0 ) then
								hydw_hydw_recip=recip_hyd_hyd_dist(test_hyd_2, test_hyd)
								term_6=term_6+v_inter_coulomb(hydw_hydw_recip,qwat_h,qwat_h)
								!!! Checked may 15 2017	
							end if
						end do  	
					end if
				end do	
				!!! ===== cl of hcl to water 3 body potential =====	
				hcl_oxw1=cl_o_dis(j)
				hcl_oxw2=cl_o_dis(k)
				hcl_3_body=hcl_3_body+rk3_hcl(hcl_oxw1,hcl_oxw2,fac3,alph2,alph2,r1sw,r2sw)
				!!! Checked July 10 2017	
			end do ! ends do loop of k=j+1 oxygens
			!!! ===== HCl and Water Potential from Wick =====	
			clhcl_oxw=cl_o_dis(j) 
			clhcl_oxw_recip=recip_cl_o_dis(j)
			hhcl_oxw_recip=recip_ox_hyd_bond(j,h_hcl)
			hhcl_oxw=ox_hyd_bond(j,h_hcl)	
			term_vrep=term_vrep+v_hcl_rep(clhcl_oxw,hhcl_oxw,repa,repb,repd,repc)  
			term_7=term_7+v_inter_lj(clhcl_oxw_recip,sig_cl_oxw_12_w,sig_cl_oxw_6_w,eps_cl_oxw_4_w)
			!!! lj cl of  hcl to oxygen of water	
			!!! Checked may 30th	
			!!! new is after
			term_7=term_7+v_inter_lj(hhcl_oxw_recip,sig_hhcl_oxw_12,sig_hhcl_oxw_6,eps_hhcl_oxw_4)
			!!! lj hydrogen of hcl to oxygen of water
			!!! Checked may 30th	
			term_8=term_8+v_inter_coulomb(clhcl_oxw_recip,qwat_o,q_clhcl)
			!!! checked june 1st	
			term_8=term_8+v_inter_coulomb(hhcl_oxw_recip,qwat_o,q_hhcl)
			!!! checked june 1st	
			do k=1,3
				test_hyd_3=wat_con(j,(3*i)-(3*hydronium_states)-3+k)
				if (test_hyd_3 .NE. 0) then 
					clhcl_hydw_recip=recip_cl_h_dis(test_hyd_3) 
					hhcl_hydw_recip=recip_hyd_hyd_dist(h_hcl,test_hyd_3)	
					term_8=term_8+v_inter_coulomb(clhcl_hydw_recip,qwat_h,q_clhcl)   
					!!! checked june 7th	
					term_8=term_8+v_inter_coulomb(hhcl_hydw_recip,qwat_h,q_hhcl)
					!!! checkedjune 7th	
					term_7=term_7+v_inter_lj(hhcl_hydw_recip,sig_hhcl_hydw_12, sig_hhcl_hydw_6,eps_hhcl_hydw_4)
					!!! lj h of hcl to hydrogen of water	
					!!! checked june 1st 2017	
					term_7=term_7+v_inter_lj(clhcl_hydw_recip,sig_clhcl_hydw_12,sig_clhcl_hydw_6,eps_clhcl_hydw_4)	
					!!! Checked june 1st 2017	
					!!! lj of cl of hcl to hydrogen of water 
				end if
			end do
		end do ! ends the do loop of j oxygens
		wat_hcl_pot=term_6+term_7+term_8+term_9+hcl_3_body+term_vrep 
		hamil_matrix(i,i)=wat_hcl_pot
		!print*, "hcl_3_body",hcl_3_body ! this is good		
		!print*, "term_6 ",term_6 ! is good 
		!print*, "term_7",term_7 ! is good	
		!print*, "term_8",term_8 ! this is good
		!print*, "term_9", term_9 ! this is good 
		!print*, "term_vrep",term_vrep !this is good  
		!print*, "total potential hcl", wat_hcl_pot ! this is good
	end if  ! end the if if its above or below hydronium states 
end do ! end loop through the number of states 
! Diagonal potential has been completely checked june 7th
!!! ===== Off-Diagonal Potential =====
!term_10 ox zundel ion to water
!term_11 shared hydrogen of zundel to water
!term_12 other zundel hydrogens to water

counter=0
do i=1,num_states-1
	if (i .LE. hydronium_states) then 
		do j=i+1,num_states
			if (j .LE. hydronium_states) then
				if (off_diag_state(i,j) .EQV. .TRUE.) then
				term_10=0.d0
        			term_11=0.d0
        			term_12=0.d0	
				temp=0.d0
					counter=counter+1
					zund_ox_1=zund_ox(counter,1)
					zund_ox_2=zund_ox(counter,2)
					zund_h_shar=zund_ox(counter,3)
					do k=1,num_ox
						if(k .NE. zund_ox_1 .AND. k .NE. zund_ox_2) then
							!!! ox zundel to ox water
							recip_zund_ox_1_oxw=recip_ox_ox_dis(zund_ox_1,k)
							recip_zund_ox_2_oxw=recip_ox_ox_dis(zund_ox_2,k)
							term_10=term_10+v_inter_coulomb(recip_zund_ox_1_oxw,qwat_o,qzund_o)
							term_10=term_10+v_inter_coulomb(recip_zund_ox_2_oxw,qwat_o,qzund_o)
							!!! Checked June 19th 2017	
							!!! shared zundel proton to ox of water
							r_sharz_hyd_oxw=recip_ox_hyd_bond(k,zund_h_shar)
							term_11=term_11+v_inter_coulomb(r_sharz_hyd_oxw,qwat_o,qzund_h_s)
							!!! Checked June 19th 2017	
							do l=1,3
								wat_h=zundel_con(k,counter*3-3+l)
								if (wat_h .NE. 0) then
									!!! zundel oxygen to water hydrogen
									r_z_ox_h_wat_1=recip_ox_hyd_bond(zund_ox_1,wat_h)
									r_z_ox_h_wat_2=recip_ox_hyd_bond(zund_ox_2,wat_h)
									term_10=term_10+v_inter_coulomb(r_z_ox_h_wat_1,qwat_h,qzund_o)
									term_10=term_10+v_inter_coulomb(r_z_ox_h_wat_2,qwat_h,qzund_o)
									!!! Checked June 19th 2017	
									!!! shared proton of zundel to hydrogen of water
									r_zsh_wath=recip_hyd_hyd_dist(zund_h_shar,wat_h)
									term_11=term_11+v_inter_coulomb(r_zsh_wath,qwat_h,qzund_h_s)
									!!! Checked June 19th 2017	
								end if
							!end do
								!!! hydrogen from zundel that's not the shared proton to the oxygen of water 
								hydron_hyd_1=zundel_con(zund_ox_1,counter*3-3+l)
								hydron_hyd_2=zundel_con(zund_ox_2,counter*3-3+l)
								if (hydron_hyd_1 .NE. zund_h_shar ) then
									rzh_oxw_1=recip_ox_hyd_bond(k,hydron_hyd_1)
									term_12=term_12+v_inter_coulomb(rzh_oxw_1,qwat_o,qzund_h)
									!!! Checked June 19th 2017	
									do m=1,3
										!!! hydrogen from zundel that's not shared to hydrogen of water
										test_hyd=zundel_con(k,counter*3-3+m)
										if (test_hyd .NE. 0) then
											rzh_hydw=recip_hyd_hyd_dist(hydron_hyd_1,test_hyd)
											term_12=term_12+v_inter_coulomb(rzh_hydw,qwat_h,qzund_h)
											!!! Checked June 19th 2017	
										end if
									end do
								end if	
								if (hydron_hyd_2 .NE. zund_h_shar) then
									rzh_oxw_2=recip_ox_hyd_bond(k,hydron_hyd_2)
									term_12=term_12+v_inter_coulomb(rzh_oxw_2,qwat_o,qzund_h)
									!!! Checked June 19 2017	
									do m=1,3
										test_hyd=zundel_con(k,counter*3-3+m)
										if (test_hyd .NE. 0) then
											rzh_hydw=recip_hyd_hyd_dist(hydron_hyd_2,test_hyd)
											term_12=term_12+v_inter_coulomb(rzh_hydw,qwat_h,qzund_h)
											!!! Checked June 19th 2017	
										end if
									end do
								end if
							end do
						end if	
					end do
					!!! A(Roo,q) part
					roo_dist=ox_ox_dis(zund_ox_1,zund_ox_2)
					q1=(xx_O(zund_ox_1)+xx_O(zund_ox_2))/2.d0-xx_H(zund_h_shar)
					q2=(yy_O(zund_ox_1)+yy_O(zund_ox_2))/2.d0-yy_H(zund_h_shar)
					q3=(zz_O(zund_ox_1)+zz_O(zund_ox_2))/2.d0-zz_H(zund_h_shar) 
					q=(q1**2)+(q2**2)+(q3**2) ! this is really q-squared
					scale_fact=arooq(roo_dist,q)
					!!! Checked June 19th 2017	
					v_ij_ex=term_10+term_11+term_12
					!!! Checked June 19th 2017	
					!print*, "v_ij_ex",v_ij_ex
					off_hyd_hyd_pot=(v_ij_const+v_ij_ex)*scale_fact
					!!! Checked June 19th 2017	
					hamil_matrix(i,j)=off_hyd_hyd_pot
					hamil_matrix(j,i)=off_hyd_hyd_pot
					!print*, "zundel state ",counter
					!print*, "temp",temp	
					!print*, "scale factor",scale_fact !this term is correct	
					!print*, "term_10",term_10 !this term is correct 
					!print*, "term_11",term_11 !this term is correct 
					!print*, "term_12", term_12 ! this term is correct 
					!print*, "off_hyd_hyd_pot",off_hyd_hyd_pot !this term is correct 
				end if
			else if (j .GT. hydronium_states) then 
				if (off_diag_state(i,j) .EQV. .TRUE.) then
					counter=counter+1
					zund_hcl_hyd=zund_ox(counter,3)
					hcl_h_zund_dist=cl_h_dis(zund_hcl_hyd)
					off_hyd_hcl_pot=hyd_hcl_zund(hcl_h_zund_dist,aof,bof,a1,b1,cons1,cons2)
					!!! Checked June 20 2017	
					!print*, "zundel state",counter
					!print*, "h ",zund_hcl_hyd	
					!print*, "rhclhyd",off_hyd_hcl_pot
					!print*, "off hcl hydronium potential ", off_hyd_hcl_pot
					hamil_matrix(i,j)=off_hyd_hcl_pot
					hamil_matrix(j,i)=off_hyd_hcl_pot
				end if
			end if
		end do
	end if
end do		

!do i=1,num_states
!	do j=1,num_states
!	print*, "at ",i," and ",j," ham is ", hamil_matrix(i,j)
!	end do
!end do
allocate(v(num_states,num_states))
allocate (d(num_states),eig_vec(num_states))
call jacobi(hamil_matrix,num_states,d,v,nrot)
! d contains the eigenvalues
! v contains the eigenvectors

!----find the lowest eigenvalue and print out the geometry
	min=1.0d10
        grstate = 3
        do i = 1,num_states
           semper = d(i)
           if (semper.lt.min) then
              min = semper
              grstate = i
              ground_state = grstate
           endif
        enddo

          !print*,'eng',min
         !stop
        do i=1,num_states
       eig_vec(i)=v(grstate,i)
	!print*, "eigenvector"
         !print*,i,eig_vec(i)
        enddo
	!print*, "ground_state"
         !print*,ground_state,d(ground_state)
           eng=d(ground_state)
            !print*,eng
          !stop
deallocate(v)
deallocate (d,eig_vec)


!!! ===== Functions =====
contains

real*8 function degree_to_radian(deg)
	implicit none
	real*8, intent(in) ::deg
	real*8 :: pi
	pi=acos(-1.d0)
	degree_to_radian=deg*(pi/180.d0)
end function

real*8 function v_hydron_intra_bond(roh)  !roh  can be given different dummy names like i and k
	implicit none
	!!! MS-EVB3 hydronium intramolecular bond terms !!!
	!!! in-term:bond distance between oxygen of hydronium and its hydrogens
	real*8, intent(in) ::roh	  
	real*8 :: d_oh, alph_oh, r_o_oh 
	!!! ===== Parameters Voth MSEVB3 ===== !!!
	!d_oh=88.96d0	 ! kcal/mol msevb 3
	!alph_oh=2.1d0	 ! inverse angstroms
	!r_o_oh=1.d0	 ! angstroms   
	!!! ===== Parameters voth MSEVB 3.2 ===== 
	d_oh=79.0864d0	! kcal/mol
	alph_oh=2.0834d0 ! inverse angstroms
	r_o_oh=0.98d0 	! angstroms
	!!! ===== Function ===== !!!
	v_hydron_intra_bond=d_oh*((1.d0-exp(-alph_oh*(roh-r_o_oh)))**2)
end function

real*8 function v_hydron_intra_ang(ang)
	implicit none
	!!! MS-EVB3 hydronium intramolecular angle terms !!!
	!!! in-term: angle between hydronium oxygen and two of its hydrogens
	real*8, intent(in) :: ang
	real*8 :: k_alph, alph_o
	!!! ===== Paramters Voth MSEVB3 ===== !!!
 	k_alph=77.4868d0 ! kcal/mol/rad_squared
        alph_o=111.7269d0 !in degrees
        alph_o=alph_o*dacos(-1.d0)/180.d0! this was converted to radians
	!!! ===== Function ===== !!!
	v_hydron_intra_ang=(1.d0/2.d0)*(k_alph*((ang-alph_o)**2))
end function

real*8 function v_inter_lj(rab_recip,sig_12, sig_6,eps_4)
	implicit none
	!!! in term: recip distance between two atoms and their sigma values raised
	!!! to the correct value and  and epsilon multiplided by 4 values 
	real*8, intent(in) :: rab_recip, sig_12, sig_6, eps_4
	!!! ===== Function ===== !!!
	v_inter_lj=(eps_4)*((sig_12*(rab_recip**12))-(sig_6*(rab_recip**6)))

end function

real*8 function v_inter_coulomb(rab, qa, qb)
	implicit none
	!!! in term: reciprocal oxygen hydrogen distance charge on particle a and particle b
	real*8, intent(in) :: rab, qa, qb
	real*8 :: coul_const
	!!! ===== Function ===== !!!
	!coul_const= 331.9208506 ! kcal*angstroms/mol*e^2
	coul_const=332.0638186797117d0
	v_inter_coulomb=coul_const*qa*qb*rab
end function

real*8 function v_inter_oo_rep_msevb3(roo,q1,q2,q3)
	implicit none
	!!! in terms: hydronium ox - wat ox dis, and thre hydronium hydrogen to wat-ox dis
	real*8, intent(in) :: roo,q1,q2,q3 
	real*8 :: b_1, b_2, b_3, doo, sw, swa, swb,swc
	!!! ===== Parameters Voth MSEVB3 =====
	!b_1=11.2600138d0	!kcal/mol
	!b_2=1.1d0		!inverse angstroms
	!b_3=2.12d0		!inverse angstroms squared
	!doo=2.4d0		!2.4 Angstroms
	!!! ===== Parameters MSEVB 3.2 =====
	b_1=9.9178410d0 	!kcal/mol
	b_2=1.1021518d0		!inverse angstroms
	b_3=2.0066249d0		! inverse angstroms squared
	doo=2.4d0
	!!! switching function values are given by msevb 3.2 
	if (roo .LT. 2.59d0) then
		sw=1.d0
	else if (roo .GE. 2.59d0 .AND. roo .LE. 2.98d0) then
		sw=1.d0-(16.85800502d0*((roo-2.59d0)**2)*(6.35d0-(2.d0*roo)))
	else if (roo .GT. 2.98d0) then
		sw=0.d0
	end if
	!!! ===== Function ===== !!!
	v_inter_oo_rep_msevb3=((b_1*exp(-b_2*(roo-doo)))*(exp(-b_3*(q1))+exp(-b_3*(q2))+exp(-b_3*(q3))))*sw
end function

real*8 function v_inter_hydo_rep_msevb3(rho)
	implicit none
	!!! in terms: the distance from the hydronium hydrogen to a water oxygen 
	real*8, intent(in) :: rho
	real*8 :: c_1, c_2, doh, sw
	!!! ===== Paramters Voth MSEVB3 =====
	!c_1=4.5715736d0	!kcal/mol
	!c_2=2.1d0	!inverse angstroms
	!doh=1.0d0	!angstroms
	!!! ===== Parameters Voth MSEVB 3.2 =====
	c_1=5.0917317d0		!kcal/mol
	c_2=8.9920023d0 	! inverse angstroms
	doh=1.0d0			! angstroms
	if (rho .LT. 1.59d0) then
		sw=1.d0
	else if (rho .GE. 1.59d0 .AND. rho .LE. 2.59d0) then
		sw=1.d0-(((rho-1.59d0)**2)*(6.18d0-(2d0*rho)))
	else if (rho .GT. 2.59d0) then
		sw=0
	end if
	!!! ===== Function ===== !!!
	v_inter_hydo_rep_msevb3=c_1*exp(-c_2*(rho-doh))*sw
end function
	
real*8 function v_intra_spcfw(roh1,roh2,ang)
	implicit none
	!!! in terms: water oh distance 1,2, and their angle
	real*8, intent(in) :: roh1,roh2,ang
	real*8 :: kb, roh_eq, ka, ang_hoh_eq
	!!! ===== Paramters SPC/FW ===== !!!
	kb=1059.162d0	!kcal/mol/angs^2
	roh_eq=1.012d0	!angstroms
	ka=75.9d0	! kcal/mol/rad^2
	ang_hoh_eq=1.97641d0	!radian
	!!! ===== Function ===== !!!
	v_intra_spcfw=((kb/2.d0)*(((roh1-roh_eq)**2)+((roh2-roh_eq)**2)))+((ka/2.d0)*((ang-ang_hoh_eq)**2))
end function 

real*8 function v_inter_bp_exp6(eps,sig,rij,recip_rij)
	implicit none
	!!! for chlorine and hcl buckingham potential
	!!! two body interactions
	!!! in terms: mixing rules epsilon, sigma, distance and recip dist
	real*8, intent(in) :: eps, sig, rij, recip_rij
	real*8 :: lam, lam_frac, sig_recip 
	!!! ===== Parameters Wick BP for ===== 	
	lam=13.5d0
	lam_frac=(6.d0)/lam  
	sig_recip=(1.d0)/sig 
	!!! ===== Function ===== !!!
	v_inter_bp_exp6=eps*((lam_frac*exp(lam*(1-(rij*sig_recip))))-((sig*recip_rij)**6)) 
end function

real*8 function v_intra_hcl_morse(r_h_cl)
	implicit none
	!!! this is the morse potential for the intramolecular hcl term
	!!! in term: distance between the h and cl
	real*8, intent(in) :: r_h_cl
	real*8 ::  d,alph,ro,uo  
	!!! ===== Parameters Wick for HCl =====
	d=106.2d0	!kcal/mol
	alph=1.913d0	!check units
	ro=1.27332d0	!angstroms
	uo=-165.11d0	! kcal/mol
	!!! ===== Function ===== !!!
	v_intra_hcl_morse=(d*((1-exp(-alph*(r_h_cl-ro)))**2)) +uo 
end function

real*8 function arooq(roo,q)
        implicit none
        !!! this is for the off-diagonal potential geometry based scalar
        !!! for the hydronium hydronium interactions
        real*8, intent(in) :: roo,q
        real*8 :: gamma,p,k,doo,beta,roo_o_1, p_p, alpha,roo_o_2
        real*8 :: term_1, term_2, term_3, term_4
        !!! ===== Parameters Voth msevb3 =====
        !gamma=1.8302895d0       ! angstroms -2
        !p=0.2327260             ! no units
        !k=9.562153d0            ! angstroms -2
        !doo=2.94d0              ! angstroms
        !beta=6.0179066d0        ! angstroms -1
        !roo_o_1=3.1d0           ! angstroms
        !p_p=10.8831327d0                ! angstroms -1
        !alpha=10.0380922d0      ! angstroms -1
        !roo_o_2=1.8136426d0     ! angstroms
	!!! ===== PARAMETERS VOTH MSEVB3.2 =====
	gamma=1.783170d0	! angstroms -2        
	p=0.1559053d0 
	k=5.0664471d0		! angstroms -2	
	doo=2.8621690d0         ! angstroms 
	beta=5.2394128d0	! angstroms -1
	roo_o_1=2.9425969d0	! angstroms
	p_p=7.6147672d0		! angstroms -1
	alpha=7.4062624d0		! angstroms -1	
	roo_o_2=1.8d0		! angstroms
	!!! ===== There are four parts to this whole function ====
        term_1= exp(-gamma*q) !!! this term isn't done yet
        term_2=1+p*exp(-k*((roo-doo)**2))
	term_3=(0.5d0)*(1-tanh(beta*(roo-roo_o_1)))
        term_4= p_p*exp(-alpha*(roo-roo_o_2))
        arooq=term_1*term_2*(term_3+term_4)
	!!! Checked June 19th 2017
end function

real*8 function hyd_hcl_zund(rh_hcl,a,b,a1,b1,cons1,cons2)
        !!! This is for the off-diagonal potential connecting states for
        !!! the hydronium hydrogen and hcl hydrogen
        !!! rh_hcl is the distance between the h of the shared proton between hydronium and hcl
        !!! and the distance is between that h and the cl
        real*8, intent(in) :: rh_hcl ,a,b,a1,b1,cons1,cons2
        real*8 :: sw, r0_s,r0_s1,sw1,ro,ro1
        !!! ===== Parameters =====
       ! a=-20.d0  ! these have to be changed
       ! b=2.5d0  ! these have to be changed
	!ro=1.4d0
!          print*,a,b
!           stop
	r0_s=2.6d0
	r0_s1=1.5d0
           ro=1.75d0
           ro1=1.4d0
        !!! ===== The Function =====
	sw= (1/2.d0)*(1-tanh(40.d0*(rh_hcl-r0_s)))   
	sw1= (1/2.d0)*(1-tanh(40.d0*(rh_hcl-r0_s1)))   
     hyd_hcl_zund=(a*exp(-b*((rh_hcl-ro)**2))+cons1)*sw+ &
&   (a1*exp(-b1*((rh_hcl-ro1)**2))+cons2)*sw1

	!print*, "sw",sw
	!print*, "sw1",sw1
	!print*, "term_1",(a*exp(-b*((rh_hcl-ro)**2))+cons1)*sw
	!print*, "term_2",(a1*exp(-b1*((rh_hcl-ro1)**2))+cons2)*sw1	
	! Checked June 20 2017
end function

real*8 function rk3(r1,r2,fac,alpha,beta)
	!!! This is the rk3 potential
	real*8, intent(in) :: alpha, fac, beta, r1, r2 
	real*8 ::sw1,sw2,rc1,rc2, g, sig, a
	!!! mw-cl potential
        rc1=3.d0
        rc2=4.d0
        !!! switching function values are given by msevb 3.2
        if (r1 .LT. rc1) then
                sw1=1.d0
        else if (r1 .GE. rc1 .AND. r1 .LE. rc2) then
                sw1=1.d0-((rc2-rc1)**(-3.d0))*((r1-rc1)**2)*(3.d0*rc2-rc1-2.d0*r1)
        else if (r1.gt.rc2) then
                sw1=0.d0
        end if
        if (r2 .LT. rc1) then
                sw2=1.d0
        else if (r2 .GE. rc1 .AND. r2 .LE. rc2) then
                sw2=1.d0-((rc2-rc1)**(-3.d0))*((r2-rc1)**2)*(3.d0*rc2-rc1-2.d0*r2)
        else if (r2.gt.rc2) then
                sw2=0.d0
        end if
	!!! ===== Function ===== !!!
             if(sw1.ne.0.d0.and.sw2.ne.0.d0)then
              rk3=fac*exp(-alpha*r1)*exp(-beta*r2)*sw1*sw2
             else
            rk3=0.d0
            endif
end function

real*8 function rk3_hyd(r1,r2,fac,alpha,beta)
        !!! This is the rk3 potential
        real*8, intent(in) :: alpha, fac, beta, r1, r2
        real*8 ::sw1,sw2,rc1,rc2,rc3,rc4, g, sig, a
        !!! mw-cl potential
        rc1=3.d0
        rc2=4.d0
	rc3=3.1
	rc4=3.5      
        if (r1 .LT. rc1) then
                sw1=1.d0
        else if (r1 .GE. rc1 .AND. r1 .LE. rc2) then
                sw1=1.d0-((rc2-rc1)**(-3.d0))*((r1-rc1)**2)*(3.d0*rc2-rc1-2.d0*r1)
        else if (r1.gt.rc2) then
                sw1=0.d0
        end if
        if (r2 .LT. rc3) then
                sw2=1.d0
        else if (r2 .GE. rc3 .AND. r2 .LE. rc4) then
                sw2=1.d0-((rc4-rc3)**(-3.d0))*((r2-rc3)**2)*(3.d0*rc4-rc3-2.d0*r2)
        else if (r2.gt.rc4) then
                sw2=0.d0
        end if
        !!! ===== Function ===== !!!
             if(sw1.ne.0.d0.and.sw2.ne.0.d0)then
              rk3_hyd=fac*exp(-alpha*r1)*exp(-beta*r2)*sw1*sw2
             else
            rk3_hyd=0.d0
            endif
end function

real*8 function v_hcl_rep(rocl,roh_hcl,repa,repb,repd,repc)
	real*8, intent(in) :: rocl,roh_hcl,repa,repb,repd,repc
	v_hcl_rep=repa*exp(-repb*(rocl-repd))*exp(-repc*(roh_hcl**2))
end function 

real*8 function v_inter_rep_hydron_cl(rclo,q1,q2,q3,b_clo_1,b_clo_2,b_clo_3,doo_clo)
        implicit none
        !!! in terms: hydronium ox - wat ox dis, and thre hydronium hydrogen to wat-ox dis
        real*8, intent(in) :: rclo,q1,q2,q3,b_clo_1,b_clo_2,b_clo_3,doo_clo 
        real*8 ::  sw 
        !!! switching function values are given by msevb 3.2
        if (rclo .LT. 2.59d0) then
                sw=1.d0
        else if (rclo .GE. 2.59d0 .AND. rclo .LE. 2.98d0) then
                sw=1.d0-(16.85800502d0*((rclo-2.59d0)**2)*(6.35d0-(2.d0*rclo)))
        else if (rclo .GT. 2.98d0) then
                sw=0.d0
        end if
        !!! ===== Function ===== !!!
        v_inter_rep_hydron_cl=((b_clo_1*exp(-b_clo_2*(rclo-doo_clo)))&
	&*(exp(-b_clo_3*(q1))+exp(-b_clo_3*(q2))+exp(-b_clo_3*(q3))))*sw
end function

real*8 function rk3_hcl(r1,r2,fac,alpha,beta,r1sw,r2sw)
        !!! This is the rk3 potential
        real*8, intent(in) :: alpha, fac, beta, r1, r2,r1sw,r2sw 
        real*8 ::sw1,sw2,rc1,rc2, g, sig, a
        !!! mw-cl potential
        rc1= r1sw 
        rc2= r2sw 
        !!! switching function values are given by msevb 3.2
        if (r1 .LT. rc1) then
                sw1=1.d0
        else if (r1 .GE. rc1 .AND. r1 .LE. rc2) then
                sw1=1.d0-((rc2-rc1)**(-3.d0))*((r1-rc1)**2)*(3.d0*rc2-rc1-2.d0*r1)
        else if (r1.gt.rc2) then
                sw1=0.d0
        end if
        if (r2 .LT. rc1) then
                sw2=1.d0
        else if (r2 .GE. rc1 .AND. r2 .LE. rc2) then
                sw2=1.d0-((rc2-rc1)**(-3.d0))*((r2-rc1)**2)*(3.d0*rc2-rc1-2.d0*r2)
        else if (r2.gt.rc2) then
                sw2=0.d0
        end if
        !!! ===== Function ===== !!!
             if(sw1.ne.0.d0.and.sw2.ne.0.d0)then
              rk3_hcl=fac*exp(-alpha*r1)*exp(-beta*r2)*sw1*sw2
             else
            rk3_hcl=0.d0
            endif
end function


end subroutine potential
