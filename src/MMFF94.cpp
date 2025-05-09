#include<iostream>
#include<string>
#include<armadillo>
#include<cmath>
#include <iomanip>
#include"math_functions.h"



/*
STRUCTS AND CLASSES USED TO REPRESENT THE MOLECULAR TOPOLOGY:
*/

struct Bond {
	double distance;
	int atom1;
	int atom2;
	// add other constants later e.g. force constant etc.
};



class Molecule {
public:
	// Attributes
	int nr_atoms;
	arma::vec atom_vec;
	arma::vec enumerated_atom_vec;
	arma::mat xyz_mat;
	std::vector<std::pair<int, arma::vec>> neighbor_vec;
	std::vector<Bond> bond_vec;
	std::vector<std::pair<double, std::array<int, 3>>> bond_angles;
	std::vector<double> improper_angles;    								 // out of plane angles
	std::vector<double> torsional_angles;
	std::vector<std::pair<double, double>> interacting_atoms; 						 // vector holding pairs that are separated by three or more bonds
	std::vector<std::pair<double, double>> atoms_1_4;							 // vector holding pairs that are separated by exactly three bonds

	// Methods
	void makeNeighborList () {
		double rH = 31 * 0.01;  					   			 // pm * 0.01 Å/pm, covalent radius H
		double rC = 76 * 0.01;  					  			 // pm * 0.01 Å/pm, covalent radius sp^3 C

		for (int i = 0; i < nr_atoms; i++) {
			std::vector<std::pair<int, arma::vec>> neighborlst_i;
			arma::vec neighbors;

			for (int j = 0; j < nr_atoms; j++) {  // two first C-H bonds, last C-C bonds
				if (euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) && 
				enumerated_atom_vec(j) != enumerated_atom_vec(i) 			        && 
				atom_vec(i) == 6 								&& 
				atom_vec(j) == 1                                     
				||				 																									    
				euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rC + rC)     && 
				enumerated_atom_vec(j) != enumerated_atom_vec(i) 				&& 
				atom_vec(i) == 6 && atom_vec(j) == 6                                           
				||	  																							            	
				euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) 	&& 
				enumerated_atom_vec(j) != enumerated_atom_vec(i) 				&& 
				atom_vec(i) == 1 								&& 
				atom_vec(j) == 6)                                 
				{
					neighbors.insert_rows(neighbors.n_rows, arma::rowvec({enumerated_atom_vec(j)}));
				}
			}
			std::pair<double, arma::vec> entry = {enumerated_atom_vec(i), neighbors};
			neighbor_vec.push_back(entry);
		}
		std::cout << "--------------------------------" << std::endl;
		std::cout << "Cartesian coordinates of atoms:" << std::endl;
		std::cout << "--------------------------------" << std::endl;
		xyz_mat.print();
		std::cout << "Number of atoms in molecule: " << nr_atoms << std::endl;
	}


	bool bondExists(double a1, double a2, const std::vector<Bond> bond_vec) {
		for (const Bond b : bond_vec) {
			if ((b.atom1 == a1 && b.atom2 == a2) || (b.atom1 == a2 && b.atom2 == a1)) {
				return true;
			}
		}
		return false;
	}



	// This method creates chemical bond structs between neighbouring atoms. TODO ADD CHECK THAT #BONDS = 3n + 1
	void makeChemicalBonds () {
		for (int i = 0; i < nr_atoms; i++) {
			for (int j : neighbor_vec[i].second) {
				if (i != j && neighbor_vec[i].second.n_elem > 1) { 					// Base the connection on the C-atoms
					int atom1 = i;
					int atom2 = j;

					if (!bondExists(atom1, atom2, bond_vec)) {
						Bond bond;
						bond.distance = euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t());
						bond.atom1 = enumerated_atom_vec(i);
						bond.atom2 = enumerated_atom_vec(j);
						bond_vec.push_back(bond);
					}
				}
			}
		}
		std::cout << "--------------------------------" << std::endl;
		std::cout << "Chemical bonds [Å]:" << std::endl;
		std::cout << "--------------------------------" << std::endl;

		for (const auto& bond : bond_vec) {
			std::cout << bond.atom1 << "-" << bond.atom2 << " " << bond.distance << std::endl;
		}

	}



	// This method computes the angles between all bonds. TODO ADD CHECK THAT #BONDANGLES = 6n
	void getBondAngles () {
		std::cout << "--------------------------------" << std::endl;
		std::cout << "Bond angles [deg]" << std::endl;
		std::cout << "--------------------------------" << std::endl;
		for (int i = 0; i < bond_vec.size(); i++) {
			for (int j = 0; j < i; j++) {
				int center_atom;
				if (bond_vec[i].atom1 == bond_vec[j].atom1) center_atom = bond_vec[i].atom1;
				else if (bond_vec[i].atom2 == bond_vec[j].atom2) center_atom = bond_vec[i].atom2;
				else if (bond_vec[i].atom1 == bond_vec[j].atom2) center_atom = bond_vec[i].atom1;
				else if (bond_vec[i].atom2 == bond_vec[j].atom1) center_atom = bond_vec[i].atom2;
				else continue;

				int external_atom1;
				int external_atom2;
				double angle;
				std::array<int, 3> atoms;

				if (center_atom == bond_vec[i].atom1) {
					external_atom1 = bond_vec[i].atom2;
				} else {
					external_atom1 = bond_vec[i].atom1;
				}

				if (center_atom == bond_vec[j].atom1) {
					external_atom2 = bond_vec[j].atom2;
				} else {
					external_atom2 = bond_vec[j].atom1;
				}

				angle = acos(arma::dot(getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom1).t()), getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom2).t()))) * (180.0/acos(-1.0));
				atoms = {external_atom1, center_atom, external_atom2};
				bond_angles.push_back(std::make_pair(angle, atoms));

				std::cout << external_atom1 << "-" << center_atom << "-" << external_atom2 << " " << angle << std::endl;
			}
		}
	}



	/*
	This method computes the OOP angles at tricoordinate centers however I did not realize that koop parameters needed to compute the energies are not given for H-C-C-C or H-C-C-H etc. becuase Out-of-plane (OOP) bending terms are only relevant for planar (spB2) centers and thus koop parameters are given for e.g. aromatics, alkenes or amides where planarity needs to be enforced. In the MMFF94 parameter set, there are no KOOP parameters for combinations like C-C-C-H or H-C-C-H involving only spB3 carbons and hydrogens. Note that is works only fo methane and ethane. It starts couble counting some combinatiosn for propane and longer due to C-C-C bonding.


		void getOOPAngles() {
	    		double improper_angle;
	    		double phijk;

	    		arma::vec e1;
	    		arma::vec e2;
	    		arma::vec e3;

	    		std::cout << "Improper angles [deg]" << std::endl;
			for (int bond = 0; bond < bond_vec.size(); bond++) {
				if (atom_vec[bond_vec[bond].atom1] == 6 && atom_vec[bond_vec[bond].atom2] == 6) {
					e1 = getUnitVector(xyz_mat.row(bond_vec[bond].atom1).t(), xyz_mat.row(bond_vec[bond].atom2).t());

					for (int neighbor1 : neighbor_vec[bond_vec[bond].atom1].second) {
	                                        e2 = getUnitVector(xyz_mat.row(bond_vec[bond].atom1).t(), xyz_mat.row(neighbor1).t());

	                                	for (int neighbor2 : neighbor_vec[bond_vec[bond].atom1].second) {

								if (neighbor1 != neighbor2 && bond_vec[bond].atom2 != neighbor1 && bond_vec[bond].atom2 != neighbor2) {
	                                                		e3 = getUnitVector(xyz_mat.row(bond_vec[bond].atom1).t(), xyz_mat.row(neighbor2).t());

					 				phijk = acos(arma::dot(e1, e2));
	                         					improper_angle = asin(arma::dot(arma::cross(e1, e2) / sin(phijk), e3));
									improper_angles.push_back(improper_angle);

									std::cout << neighbor1 << "-" << bond_vec[bond].atom1 << "-" << bond_vec[bond].atom2 << "-" << neighbor2 << " " << improper_angle * (180.0/acos(-1.0)) << std::endl;
								}
	                                               	}
	                                        }
					}

					e1 = getUnitVector(xyz_mat.row(bond_vec[bond].atom2).t(), xyz_mat.row(bond_vec[bond].atom1).t());

					for (int neighbor1 : neighbor_vec[bond_vec[bond].atom2].second) {
	                                        e2 = getUnitVector(xyz_mat.row(bond_vec[bond].atom2).t(), xyz_mat.row(neighbor1).t());

	                                        for (int neighbor2 : neighbor_vec[bond_vec[bond].atom2].second) {
	                                                if (neighbor1 != neighbor2 && bond_vec[bond].atom1 != neighbor1 && bond_vec[bond].atom1 != neighbor2) {
	                                                        e3 = getUnitVector(xyz_mat.row(bond_vec[bond].atom2).t(), xyz_mat.row(neighbor2).t());

	                                                        phijk = acos(arma::dot(e1, e2));
	                                                        improper_angle = asin(arma::dot(arma::cross(e1, e2) / sin(phijk), e3));
	                                                        improper_angles.push_back(improper_angle);

	                                        std::cout << neighbor1 << "-" << bond_vec[bond].atom2 << "-" << bond_vec[bond].atom1 << "-" << neighbor2 << " " << improper_angle * (180.0/acos(-1.0)) << std::endl;
	                                                }
	                                        }
					}
				}
			}
	*/


bool isBonded(double a1, double a2) {
       for (const Bond b : bond_vec) {
               if ((a1 == b.atom1 && a2 == b.atom2) || (a1 == b.atom2 && a2 == b.atom1)) {
                       return true;
               }
       }
       return false;
       }

void getTorsionalAngles () {
		double torsional_angle; 
		double bond_angle_tol = 1e-8;
		double phijk;
		double phjkl;

		arma::vec e1; 
		arma::vec e2;
		arma::vec e3; 
		arma::vec cross_e1e2;
		arma::vec cross_e2e3;
	
		std::cout << "--------------------------------" << std::endl;
		std::cout << "Torsional angles [deg]" << std::endl;
		std::cout << "--------------------------------" << std::endl;

		for (int bond = 0; bond < bond_vec.size(); bond++) {												// Identify C-C bond
			if (atom_vec[bond_vec[bond].atom1] == 6 && atom_vec[bond_vec[bond].atom2] == 6) {
				for (int neighbor1 : neighbor_vec[bond_vec[bond].atom1].second) {
					for (int neighbor2 : neighbor_vec[bond_vec[bond].atom2].second) {
						if (atom_vec[neighbor1] != 6 && atom_vec[neighbor2] != 6) {   							// Identify H-C-C-H structure 
							e1 = getUnitVector(xyz_mat.row(neighbor1).t(), xyz_mat.row(bond_vec[bond].atom1).t());
							e2 = getUnitVector(xyz_mat.row(bond_vec[bond].atom1).t(), xyz_mat.row(bond_vec[bond].atom2).t());
							e3 = getUnitVector(xyz_mat.row(bond_vec[bond].atom2).t(), xyz_mat.row(neighbor2).t());


							phijk = acos(arma::dot(-e1, e2));
							phjkl = acos(arma::dot(-e2, e3));

							cross_e1e2 = arma::cross(e1, e2);
							cross_e2e3 = arma::cross(e2, e3);

							double argument = arma::dot(cross_e1e2, cross_e2e3) * std::pow(sin(phijk) * sin(phjkl), -1);			
							argument = std::max(-1.0, std::min(1.0, argument)); 						// Had to clamp the argument because nan values appearing instead of 180deg due to the argument being -1.00000000001																			

							torsional_angle = acos(argument); 
	
		
							if (sin(phijk) >= bond_angle_tol or sin(phjkl) >= bond_angle_tol) { 					// Skip unphysical torsions	
								torsional_angles.push_back(torsional_angle);							// Append in rad because will take the cosine of these later. 
								std::cout << neighbor1 << "-" << bond_vec[bond].atom1 << "-" << bond_vec[bond].atom2 << "-" << neighbor2 << " " << torsional_angle * (180.0/acos(-1.0)) << std::endl;
							}	
						} 

						if (atom_vec[neighbor1] == 6 && atom_vec[neighbor2] == 6 && !isBonded(neighbor1, neighbor2)) {							// Identify a C-C-C-C structure 
							e1 = getUnitVector(xyz_mat.row(neighbor1).t(), xyz_mat.row(bond_vec[bond].atom1).t());
                                                        e2 = getUnitVector(xyz_mat.row(bond_vec[bond].atom1).t(), xyz_mat.row(bond_vec[bond].atom2).t());
                                                        e3 = getUnitVector(xyz_mat.row(bond_vec[bond].atom2).t(), xyz_mat.row(neighbor2).t());


                                                        phijk = acos(arma::dot(-e1, e2));
                                                        phjkl = acos(arma::dot(-e2, e3));

                                                        cross_e1e2 = arma::cross(e1, e2);
                                                        cross_e2e3 = arma::cross(e2, e3);

                                                        double argument = arma::dot(cross_e1e2, cross_e2e3) * std::pow(sin(phijk) * sin(phjkl), -1);
                                                        argument = std::max(-1.0, std::min(1.0, argument));                                             // Had to clamp the argument because nan values appearing instead of 180deg due to the argument being -1.00000000001

                                                        torsional_angle = acos(argument);


                                                        if (sin(phijk) >= bond_angle_tol or sin(phjkl) >= bond_angle_tol) {                                     // Skip unphysical torsions
                                                                torsional_angles.push_back(torsional_angle);                                                    // Append in rad because will take the cosine of these later.
                                                                std::cout << neighbor1 << "-" << bond_vec[bond].atom1 << "-" << bond_vec[bond].atom2 << "-" << neighbor2 << " " << torsional_angle * (180.0/acos(-1.0)) << std::endl;
                                                        }

						}

					}
				}
			}
		}
	}



// This function computes the pairs of 1-4 interacting atoms and beyond i.e. pairs of atoms with three or more bonds between them for which vdw and electrostatic interactions become relevant.
	void getInteractingAtoms() {
		std::pair<double, double> interacting_pair;
		std::pair<double, double> one_four_interacting_pair;
		arma::vec neighbors_of_neighbors;												// Vector containing the neighbors of neighbors of an atom.

		std::cout << "--------------------------------" << std::endl;
		std::cout << "1-4 interacting pairs" << std::endl;
		std::cout << "--------------------------------" << std::endl;

		for (int i = 0; i < atom_vec.size(); i++) {											// Loop over all unique pairs
			for (int j = i + 1; j < atom_vec.size(); j++) {
				if (!arma::any(neighbor_vec[i].second == j)) {									// Skip neighbors of neighbors
					neighbors_of_neighbors.reset();

					for (int l : neighbor_vec[i].second) {
						neighbors_of_neighbors = arma::join_vert(neighbors_of_neighbors, neighbor_vec[l].second);       // Detemine all neighbors of neighbors and exclude them too.
					}

					if (!arma::any(neighbors_of_neighbors == static_cast<double>(j))) {
						interacting_pair = {i, j};
						std::cout << i << " " << j << std::endl;
						interacting_atoms.push_back(interacting_pair);
					}
				}
			}
		}
	std::cout << "--------------------------------" << std::endl;
	}
};



/*
FUNCTIONS USED TO CONSTRUCT THE FORCEFIELD AND COMPUTE THE ENERGY:
*/

// Function used to read in atomic positions from a .txt file to an armadillo matrix.
arma::mat readFile(std::string filename) {
	std::ifstream file (filename);
	if (file.is_open()) {
		std::string first_line;
		std::getline(file, first_line);  // read and discard the first line

		int nr_atoms = std::stoi(first_line);
		arma::mat input_data_matrix = arma::mat(nr_atoms, 4);

		for (int i = 0; i < nr_atoms; i++) {
			for (int j = 0; j < 4; j++) {
				file >> input_data_matrix(i,j);
			}
		}
		file.close();
		return input_data_matrix;
	} else {
		std::cerr << "Error: Failed to open file." << std::endl;
	}
}





// Function used to construct the molecular topology
Molecule buildMolecule(arma::mat input_data) {
	Molecule molecule;
	molecule.nr_atoms = input_data.n_rows;
	molecule.atom_vec = input_data.col(0);
	molecule.enumerated_atom_vec = enumerateVector(input_data.col(0));
	molecule.xyz_mat  = input_data.cols(1,3);
	molecule.makeNeighborList();
	molecule.makeChemicalBonds();
	molecule.getBondAngles();
	//molecule.getOOPAngles();
	molecule.getTorsionalAngles();
	molecule.getInteractingAtoms();
	return molecule;
}



// Function used to evaluate the bond stretching energy
double bondStretchingEnergy(Bond bond, Molecule mol) {
	double energy;
	double kb;
	double distance;
	double refdistance;
	double del_r;
	double cs; 																				// Cubic stretch constant

	if (mol.atom_vec[bond.atom1] == 6 && mol.atom_vec[bond.atom2] == 6) {
		kb = 4.258;
		refdistance = 1.508;
	}																					// C-C

	if ((mol.atom_vec[bond.atom1] == 6 && mol.atom_vec[bond.atom2] == 1) || (mol.atom_vec[bond.atom1] == 1 && mol.atom_vec[bond.atom2] == 6)) {				// C-H
		kb = 4.766;
		refdistance = 1.093;
	}

	del_r = bond.distance - refdistance;
	cs = 2;

	energy = 143.9325 * kb * 0.5 * std::pow(del_r, 2) * (1 + cs * del_r + 7 * std::pow(12, -1) * std::pow(cs, 2) * std::pow(del_r, 2));

	return energy;
}



// Function used to evaluate the angle bending energy
double angleBendingEnergy(double angle, std::array<int, 3> atoms, Molecule mol) {
	double I = atoms[0];
	double J = atoms[1];
	double K = atoms[2];
	double theta = angle;
	double angle_tol = 1;

	double energy;
	double ka;
	double reftheta;
	double del_theta;
	double cb;

	if ((mol.atom_vec[I] == 6 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 1) || (mol.atom_vec[I] == 1 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 6)) { // C-C-H or H-C-C
		ka = 0.636;
		reftheta = 110.549;
	}

	if (mol.atom_vec[I] == 1 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 1) {									       	// H-C-H
		ka = 0.516;
		reftheta = 108.836;
	}

	if (mol.atom_vec[I] == 6 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 6) {                                                                        	// C-C-C
		ka = 0.851;
		reftheta =  109.608;
	}

	del_theta = theta - reftheta;
	cb = -0.007;

	if (std::abs(theta - 180.0) <= angle_tol) {														// MMFF treats linear and non linear bonds using a different equations for the angle bending energy.
		energy = 143.9325 * ka * (1 + cos(theta));
		std::cout << "we're doomed " << std::endl;
	} else {
		energy = 0.043844 * ka * 0.5 * std::pow(del_theta, 2) * (1 + cb * del_theta);
	}

	return energy;
}



// Function used to evaluate the stretch-bend energy
double stretchBendEnergy(double angle, std::array<int, 3> atoms, Bond bond, Molecule mol) {
	double energy;
	double kIJK;
	double kKJI;
	double rij;
	double rkj;
	double refdistanceij;
	double refdistancekj;
	double del_rij;
	double del_rkj;
	double theta;
	double reftheta;
	double del_theta;

	if (mol.atom_vec[atoms[0]] == 6 && mol.atom_vec[atoms[1]] == 6 && mol.atom_vec[atoms[2]] == 6) {														// C-C-C
		kIJK = 0.206;
		kKJI = 0.206;
		refdistanceij = 1.508;
		refdistancekj = 1.508;
		reftheta = 109.608;
	}

	if (mol.atom_vec[atoms[0]] == 6 && mol.atom_vec[atoms[1]] == 6 && mol.atom_vec[atoms[2]] == 1) {               												 	// C-C-H
		kIJK = 0.227;
		kKJI = 0.070;
		refdistanceij = 1.508;
		refdistancekj = 1.093;
		reftheta = 110.549;
	}

	if (mol.atom_vec[atoms[0]] == 1 && mol.atom_vec[atoms[1]] == 6 && mol.atom_vec[atoms[2]] == 6) {                                                                                                                // H-C-C
		kIJK = 0.070;
		kKJI = 0.227;
		refdistanceij = 1.093;
		refdistancekj = 1.508;
		reftheta = 110.549;

	}

	if (mol.atom_vec[atoms[0]] == 1 && mol.atom_vec[atoms[1]] == 6 && mol.atom_vec[atoms[2]] == 1) {                                        									// H-C-H
		kIJK = 0.115;
		kKJI = 0.115;
		refdistanceij = 1.093;
		refdistancekj = 1.093;
		reftheta = 108.836;
	}

	rij = euclideanDistance(mol.xyz_mat.row(atoms[1]).t(), mol.xyz_mat.row(atoms[0]).t());
	rkj = euclideanDistance(mol.xyz_mat.row(atoms[1]).t(), mol.xyz_mat.row(atoms[2]).t());


	for (int i = 0; i < mol.bond_angles.size(); i++) {
		if (mol.bond_angles[i].second[0] == atoms[0] && mol.bond_angles[i].second[1] == atoms[1] && mol.bond_angles[i].second[2] == atoms[2]) {
			theta = mol.bond_angles[i].first;
		}
	}

	del_rij   = rij - refdistanceij;
	del_rkj   = rkj - refdistancekj;
	del_theta = theta - reftheta;

	energy = 2.51210 * (kIJK * del_rij + kKJI * del_rkj) * del_theta;

	return energy;
}



// Function used to evaluate the torsion interaction energy
double torsionEnergy(double torsional_angle, Molecule mol) {
	double energy;
	double reference_point = 0.0; 													// Reference such that the torsions be positive
	double V1 = 0.284;															// Since only linear alkanes are being considered H-C-C-H is the only possible combination that can form a torsional angle.
	double V2 = -1.386;
	double V3 = 0.314;

	energy = 0.5 * (V1 * (1 + cos(torsional_angle)) + V2 * (1 - cos(2 * torsional_angle)) + V3 * (1 + cos(3 * torsional_angle))) + reference_point;

	return energy;
}



// Function used to evaluate the Van der Waals energy
double vdwEnergy(std::pair<double, double> atoms, Molecule mol) {
	double energy;
	double Rij;
	double gammaIJ;
	double epsilonIJ;
	double RIJstar;
	double RIIstar;
	double RJJstar;
	double alphaI;
	double alphaJ;
	double AI;
	double AJ;
	double NI;
	double NJ;
	double GI;
	double GJ;

	if (mol.atom_vec[atoms.first] == 6 && mol.atom_vec[atoms.second] == 6) {   			// C - - C interaction
		alphaI = 1.050;
		alphaJ = 1.050;
		AI = 3.890;
		AJ = 3.890;
		NI = 2.490;
		NJ = 2.490;
		GI = 1.282;
		GJ = 1.282;
	}

	if ((mol.atom_vec[atoms.first] == 6 && mol.atom_vec[atoms.second] == 1) || (mol.atom_vec[atoms.first] == 1 && mol.atom_vec[atoms.second] == 6)) {                        // C - - H interaction
		alphaI = 1.050;
		alphaJ = 0.250;
		AI = 3.890;
		AJ = 4.200;
		NI = 2.490;
		NJ = 0.800;
		GI = 1.282;
		GJ = 1.209;
	}

	if (mol.atom_vec[atoms.first] == 1 && mol.atom_vec[atoms.second] == 1) {                        // H - - H interaction
		alphaI = 0.250;
		alphaJ = 0.250;
		AI = 4.200;
		AJ = 4.200;
		NI = 0.800;
		NJ = 0.800;
		GI = 1.209;
		GJ = 1.209;
	}


	RIIstar = AI * std::pow(alphaI, 0.25);
	RJJstar = AJ * std::pow(alphaJ, 0.25);
	gammaIJ = (RIIstar - RJJstar) * std::pow((RIIstar + RJJstar), -1);
	RIJstar = 0.5 * (RIIstar + RJJstar) * (1 + 0.2 * (1 - std::exp(-12 * std::pow(gammaIJ, 2))));
	epsilonIJ = 181.16 * GI * GJ * alphaI * alphaJ * std::pow(std::pow(alphaI * std::pow(NI, -1), 0.5) + std::pow(alphaJ * std::pow(NJ, -1), 0.5), -1) * std::pow(RIJstar, -6);

	Rij = euclideanDistance(mol.xyz_mat.row(atoms.first).t(), mol.xyz_mat.row(atoms.second).t());
	energy = epsilonIJ * std::pow(1.07 * RIJstar * std::pow((Rij + 0.07 * RIJstar), -1), 7) * (1.12 * std::pow(RIJstar, 7) * std::pow((std::pow(Rij, 7) + 0.12 * std::pow(RIJstar, 7)), -1) - 2);

	return energy;
}


// Function used to evaluate the energy stemming from electrostatic repulsions.
double electrostaticEnergy(Molecule mol) {
	double energy = 0.0;
	double delta = 0.05;
	int n = 1;
	double D = 1;
	double Rij;

	std::vector<double> charges(mol.atom_vec.size(), 0.0);

	// First calculate charges for all atoms
	for (int i = 0; i < mol.atom_vec.size(); i++) {
		double qi = 0.0;
		double rhoi;

		if (mol.atom_vec[i] == 6) {
			rhoi = 0.0;      // C
		} else if (mol.atom_vec[i] == 1) {
			rhoi = -0.023;   // H
		} else {
			rhoi = 0.0;      // Default for other atom types
		}

		for (int neighbor : mol.neighbor_vec[i].second) {
			double rhoK;
			if (mol.atom_vec[neighbor] == 6) {
				rhoK = 0.0;      // C
			} else if (mol.atom_vec[neighbor] == 1) {
				rhoK = -0.023;   // H
			} else {
				std::cout << "Error unknwon atomtype encountered when assigning partial charges" << std::endl;
			}
			qi += (rhoi - rhoK);
		}

		charges[i] = qi;
	}

	// Now calculate energy for the specified interacting atom pairs
	for (int i = 0; i < mol.interacting_atoms.size(); i++) {
		int atom1 = mol.interacting_atoms[i].first;
		int atom2 = mol.interacting_atoms[i].second;

		Rij = euclideanDistance(mol.xyz_mat.row(atom1).t(), mol.xyz_mat.row(atom2).t());
		energy += 332.0716 * charges[atom1] * charges[atom2] / (D * std::pow((Rij + delta), n));
	}

	return energy;
}



// Main function that executes the program.
int main(int argc, char** argv) {
	// Set decimals:
	std::cout << std::fixed << std::setprecision(7);


	// Read data:
	std::string filepath;

	if (argc < 2) {
		std::cerr << "Error: Format for input filepath must be:  " << argv[0] << " <file_path>\n";
		return 1;
	} else {
		filepath = argv[1];
	}

	arma::mat input_mat = readFile(filepath);

	// Build molecule
	Molecule mol;
	mol = buildMolecule(input_mat);


	// Compute bond energy:
	double bond_stretching_energy;
	for (int i = 0; i < mol.bond_vec.size(); i++) {
		bond_stretching_energy += bondStretchingEnergy(mol.bond_vec[i], mol);
	}


	// Compute angle bending energy:
	double angle_bending_energy;
	for (int i = 0; i < mol.bond_angles.size(); i++) {
		angle_bending_energy += angleBendingEnergy(mol.bond_angles[i].first, mol.bond_angles[i].second, mol);
	}


	// Compute stretch-bend interactions:
	double stretch_bending_energy;
	for (int i = 0; i < mol.bond_angles.size(); i++) {
		stretch_bending_energy += stretchBendEnergy(mol.bond_angles[i].first, mol.bond_angles[i].second, mol.bond_vec[i], mol);
	}

	// Compute torsional energy:
	double torsional_energy;
	for (int i = 0; i < mol.torsional_angles.size(); i++) {
		torsional_energy += torsionEnergy(mol.torsional_angles[i], mol);
	}

	// Compute the VdWEnergy:
	double VdW_energy;
	for (int i = 0; i < mol.interacting_atoms.size(); i++) {
		VdW_energy += vdwEnergy(mol.interacting_atoms[i], mol);
	}

	// Compute electrostatic energy:
	double electrostatic_energy;
	electrostatic_energy = electrostaticEnergy(mol);


	double mmff94energy; 
	mmff94energy = bond_stretching_energy + angle_bending_energy + stretch_bending_energy + torsional_energy + VdW_energy; // + electrostatic_energy; 
	

	// Output:	
	std::cout << std::endl;   
	std::cout << std::endl;   
	std::cout << std::endl;   
	std::cout << std::string(50, '-') << std::endl;
   	std::cout << "Molecular Information:" << std::endl;
   	std::cout << std::string(50, '-') << std::endl;
	std::cout << std::left  << std::setw(35) << "Number of bonds:"
                  << std::right << std::setw(10) << mol.bond_vec.size() << std::endl;
    	std::cout << std::left  << std::setw(35) << "Number of angles:"
              	  << std::right << std::setw(10) << mol.bond_angles.size() << std::endl;
    	std::cout << std::left  << std::setw(35) << "Number of interacting pairs:"
                  << std::right << std::setw(10) << mol.interacting_atoms.size() << std::endl;
    	std::cout << std::string(50, '-') << std::endl;
	std::cout << std::endl;   
	std::cout << std::endl;   
	std::cout << std::endl;   
	std::cout << std::string(50, '-') << std::endl;
    	std::cout << "Energetics:" << std::endl;
    	std::cout << std::string(50, '-') << std::endl;
    	std::cout << std::left  << std::setw(35) << "Term"
                  << std::right << std::setw(10) << "Energy" << std::endl;
    	std::cout << std::string(50, '-') << std::endl;
    	std::cout << std::left  << std::setw(35) << "bond_stretching_energy [kcal/mol]"
                  << std::right << std::setw(10) << bond_stretching_energy << std::endl;
    	std::cout << std::left  << std::setw(35) << "angle_bending_energy [kcal/mol]"
                  << std::right << std::setw(10) << angle_bending_energy << std::endl;
    	std::cout << std::left  << std::setw(35) << "stretch_bending_energy [kcal/mol]"
                  << std::right << std::setw(10) << stretch_bending_energy << std::endl;
    	std::cout << std::left  << std::setw(35) << "torsional_energy [kcal/mol]"
                  << std::right << std::setw(10) << torsional_energy << std::endl;
    	std::cout << std::left  << std::setw(35) << "vdw_energy [kcal/mol]"
                  << std::right << std::setw(10) << VdW_energy << std::endl;
    	std::cout << std::left  << std::setw(35) << "electrostatic_energy [kcal/mol]"
                  << std::right << std::setw(10) << electrostatic_energy << std::endl;
    	std::cout << std::string(50, '-') << std::endl;
	std::cout << std::left  << std::setw(35) << "MMFF94 Energy [kcal/mol]"
                  << std::right << std::setw(10) << mmff94energy << std::endl;

	return 0;
}


/*
        // Torsion debugging:
        std::cout << "TORSION DEBUGGING" << std::endl;
        const double step = 0.1;
                for (double theta = 0.0; theta <= 2 * M_PI; theta += step) {
                double energy = torsionEnergy(theta, mol);
                std::cout << theta << "," << energy << "\n";
                }

        // DEBUGGING:
        std::cout << "Atom connectivity (with element types):" << std::endl;
        for (int i = 0; i < mol.nr_atoms; ++i) {
                std::string atom_type;
                if (mol.atom_vec(i) == 6) atom_type = "C";
                else if (mol.atom_vec(i) == 1) atom_type = "H";
                else atom_type = "X"; // Unknown atom

                std::cout << "Atom " << i << " (" << atom_type << ") is bonded to: ";

                std::vector<int> connected_atoms;
                for (const Bond& bond : mol.bond_vec) {
                        if (bond.atom1 == i) {
                                connected_atoms.push_back(bond.atom2);
                        } else if (bond.atom2 == i) {
                                connected_atoms.push_back(bond.atom1);
                        }
                }

                if (connected_atoms.empty()) {
                        std::cout << "none";
                } else {
                        for (size_t j = 0; j < connected_atoms.size(); ++j) {
                                int neighbor_idx = connected_atoms[j];
                                std::string neighbor_type;
                                if (mol.atom_vec(neighbor_idx) == 6) neighbor_type = "C";
                                else if (mol.atom_vec(neighbor_idx) == 1) neighbor_type = "H";
                                else neighbor_type = "X";

                                std::cout << neighbor_idx << " (" << neighbor_type << ")";
                                if (j < connected_atoms.size() - 1) std::cout << ", ";
                        }
                }
                std::cout << std::endl;
        }
*/


