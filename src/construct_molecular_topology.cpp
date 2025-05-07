#include<iostream> 
#include<string> 
#include<armadillo>
#include<cmath>
#include<set>
#include "math_functions.h"

/*
STRUCTS AND CLASSES USED TO CONSTRUCT THE MOLECULAR TOPOLOGY:
*/
struct Bond {
	double distance; 
	double atom1;
	double atom2; 
	// add other constants later e.g. force constant etc. 
};


class Molecule {
public:
        // Attributes
        int nr_atoms;
        arma::vec atom_vec;
	arma::vec enumerated_atom_vec; 
        arma::mat xyz_mat;
	std::vector<std::pair<double, arma::vec>> neighbor_vec;
	std::vector<Bond> bond_vec; 
	std::vector<std::pair<double, arma::vec>> bond_angles; 
	std::vector<double> improper_angles;    								 // out of plane angles 
	std::vector<double> torsional_angles; 
        std::vector<std::pair<double, double>> interacting_atoms; 

        // Methods
        void makeNeighborList () {
                double rH = 31 * 0.01;  					   			 // pm * 0.01 Å/pm, covalent radius H
                double rC = 76 * 0.01;  					  			 // pm * 0.01 Å/pm, covalent radius sp^3 C

                for (int i = 0; i < nr_atoms; i++) {
                        std::vector<std::pair<double, arma::vec>> neighborlst_i;
                        arma::vec neighbors; 

                        for (int j = 0; j < nr_atoms; j++) {  // two first C-H bonds, last C-C bonds
                                if (euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 6 && atom_vec(j) == 1                                     || 																									    euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rC + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 6 && atom_vec(j) == 6                                           ||  																							            euclideanDistance(xyz_mat.row(i).t(), xyz_mat.row(j).t()) < 1.3 * (rH + rC) && enumerated_atom_vec(j) != enumerated_atom_vec(i) && atom_vec(i) == 1 && atom_vec(j) == 6)                                 {                       			
				    neighbors.insert_rows(neighbors.n_rows, arma::rowvec({enumerated_atom_vec(j)}));  
                                 }
                        }
                        std::pair<double, arma::vec> entry = {enumerated_atom_vec(i), neighbors};
                        neighbor_vec.push_back(entry);
                }
	std::cout << "Number of atoms in molecule: " << nr_atoms << std::endl;
	std::cout << "Cartesian coordinates of atoms:" << std::endl;
	xyz_mat.print(); 
        }


	bool bondExists(double a1, double a2, const std::vector<Bond> bond_vec) {	 
        for (const Bond b : bond_vec) {
                if ((b.atom1 == a1 && b.atom2 == a2) || (b.atom1 == a2 && b.atom2 == a1)) {
                        return true;
                }
        }
        return false;
	}	

	/*
	bool isBonded(double a1, double a2) {
	for (const Bond b : bond_vec) {
		if ((a1 == b.atom1 && a2 == b.atom2) || (a1 == b.atom2 && a2 == b.atom1)) {
			return true;
		}
	}
	return false; 
	}
	*/

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
	std::cout << "Chemical bonds [Å]:" << std::endl;
	
	for (const auto& bond : bond_vec) {
    		std::cout << bond.atom1 << "-" << bond.atom2 << " " << bond.distance << std::endl;
		}

	}


	void getBondAngles () {
		std::cout << "Bond angles [deg]" << std::endl;
		for (int i = 0; i < bond_vec.size(); i++) {
			for (int j = 0; j < i; j++) {
				double center_atom; 
				if (bond_vec[i].atom1 == bond_vec[j].atom1) center_atom = bond_vec[i].atom1; 
				else if (bond_vec[i].atom2 == bond_vec[j].atom2) center_atom = bond_vec[i].atom2;
				else if (bond_vec[i].atom1 == bond_vec[j].atom2) center_atom = bond_vec[i].atom1;
				else if (bond_vec[i].atom2 == bond_vec[j].atom1) center_atom = bond_vec[i].atom2;
				else continue;  
					
					double external_atom1;
					double external_atom2;
					double angle;
					arma::vec atoms; 
		
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
								
					std::cout << external_atom1 << "-" << center_atom << "-" << external_atom2 << " " << angle * (180.0/acos(-1.0)) << std::endl;
				}
			}
		}


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
	
		std::cout << "Torsional angles [deg]" << std::endl;
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
					}
				}
			}
		}
	}



	void getInteractingAtoms() {
                std::pair<double, double> interacting_pair;

		arma::vec test2;

                std::cout << "1-4 interacting pairs" << std::endl;
                for (int i = 0; i < atom_vec.size(); i++) {
                        for (int j = i+1; j < atom_vec.size(); j++) {                                                                                   // Loop over all unique pairs
                                if (!arma::any(neighbor_vec[i].second == j)) {
					for (int l : neighbor_vec[i].second) { 
                                        	test2 = arma::join_vert(test2, neighbor_vec[l].second);                					// Exclude neighbors and neighbors of neighbors
                                        }
						
                                        for (int k : neighbor_vec[i].second) {
                                                if (atom_vec[k] == 6 && !arma::any(neighbor_vec[k].second == j) && !arma::any(test2 == static_cast<double>(j))) {
                                                        interacting_pair = {i, j};
                                                        std::cout << i << " " << j << std::endl;
                                                        interacting_atoms.push_back(interacting_pair);
                                                }
                                        }
					test2.reset();
                                }
                        }
                }
	std::cout << "Number of interacting pairs: " << interacting_atoms.size() << std::endl;
        }
};


/*
FUNCTIONS USED TO CONSTRUCT THE MOLECULAR TOPOLOGY:
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


// Function used to enumerate an arma::vec object. 
arma::vec enumerateVector(arma::vec vector) {
	for (int i = 0; i < vector.n_elem; i++) {
		vector(i) = i;
	}
	return vector; 
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
        molecule.getOOPAngles();
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
double angleBendingEnergy(double angle, arma::vec atoms, Molecule mol) {
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

	if ((mol.atom_vec[I] == 6 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 1) || (mol.atom_vec[I] == 1 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 6)) {        // C-C-H or H-C-C
		ka = 0.636; 
		reftheta = 110.549; 
	}

	if (mol.atom_vec[I] == 1 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 1) {									       // H-C-H
		ka = 0.516;
		reftheta = 108.836; 
        }

	if (mol.atom_vec[I] == 6 && mol.atom_vec[J] == 6 && mol.atom_vec[K] == 6) {                                                                        // C-C-C
		ka = 0.851; 
		reftheta =  109.608; 
        } 

	del_theta = theta - reftheta; 
	cb = -0.007; 

	if (std::abs(theta - 180.0) <= angle_tol) {												// MMFF treats linear and non linear bonds using a different equations for the angle bending energy. 
		energy = 143.9325 * ka * (1 + cos(theta));
	} else {
		energy = 0.043844 * ka * 0.5 * std::pow(del_theta, 2) * (1 + cb * del_theta); 
	}

	return energy;
}






// Function used to evaluate the strecth-bend energy 




// Function used to evlaute the out of plane bending energy at triicoordinate centers 




// Function used to evaluate the torsion interaction energy
double torsionEnergy(double torsional_angle, Molecule mol) {
	double energy; 
	double V1 = 0.284;															// Since only linear alkanes are being considered H-C-C-H is the only possible combination that can form a torsional angle. 
	double V2 = -1.386;
	double V3 = 0.314;

	energy = 0.5 * (V1 * (1 + cos(torsional_angle)) + V2 * (1 - cos(2 * torsional_angle)) + V3 * (1 + cos(3 * torsional_angle)));
	
	return energy; 
}



// Function used to evaluate the van der waals energy steming from van der Waals interactions




// Function used to evaluate the energy steming from electrostatic repulsions. 



// Main function that executes the program. 
int main(int argc, char** argv) {
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

	std::cout << "In main" << std::endl; 
	std::cout << "mol.nr_atoms" << mol.nr_atoms << std::endl;
	std::cout << "mol.atom_vec" << std::endl;
	mol.atom_vec.print();
	std::cout << "mol.enumerated_atom_vec" <<  std::endl;
	mol.enumerated_atom_vec.print(); 
	std::cout << "xyz_mat" << std::endl;
	mol.xyz_mat.print();
	std::cout << "mol.neighbor_vec.size()" << mol.neighbor_vec.size() << std::endl; 
	std::cout << "mol.bond_vec.size()" << mol.bond_vec.size() << std::endl;
	std::cout << "mol.bond_angles.size()"<< mol.bond_angles.size() << std::endl;
	std::cout << "mol.improper_angles.size()" << mol.improper_angles.size() << std::endl;
	std::cout << "mol.torsional_angles.size()" << mol.torsional_angles.size() << std::endl;
	std::cout << "mol.interacting_atoms.size()" << mol.interacting_atoms.size() << std::endl;


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


	// Compute torsional energy:
	double torsional_energy;
	for (int i = 0; i < mol.torsional_angles.size(); i++) {
		torsional_energy += torsionEnergy(mol.torsional_angles[i], mol); 
	}


	std::cout << "bond_stretching_energy  [kJ/mol]" << " " << bond_stretching_energy << std::endl; 
	std::cout << "angle_bending_energy    [kJ/mol]" << " " << angle_bending_energy   << std::endl;
	std::cout << "torsional_energy        [kJ/mol]" << " " << torsional_energy	 << std::endl;


return 0;
}
