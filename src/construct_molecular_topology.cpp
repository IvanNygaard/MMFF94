#include<iostream> 
#include<string> 
#include<armadillo>
#include<cmath>
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
	std::vector<double> bond_angles; 
	std::vector<double> wilson_angles;    								 // out of plane angles 
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


	void makeChemicalBonds () {   // Want to loop over all the carbons and create bond from each carbon to its neighbour (but not doulbe count C-C bonds!)
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
				int center_atom; 
				if (bond_vec[i].atom1 == bond_vec[j].atom1) center_atom = bond_vec[i].atom1; 
				else if (bond_vec[i].atom2 == bond_vec[j].atom2) center_atom = bond_vec[i].atom2;
				else if (bond_vec[i].atom1 == bond_vec[j].atom2) center_atom = bond_vec[i].atom1;
				else if (bond_vec[i].atom2 == bond_vec[j].atom1) center_atom = bond_vec[i].atom2;
				else continue;  
					
					double external_atom1;
					double external_atom2;
					double angle;
		
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

					angle = acos(arma::dot(getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom1).t()), getUnitVector(xyz_mat.row(center_atom).t(), xyz_mat.row(external_atom2).t()))); 
	                		bond_angles.push_back(angle);  
								
					std::cout << external_atom1 << "-" << center_atom << "-" << external_atom2 << " " << angle * (180.0/acos(-1.0)) << std::endl;
				}
			}
		}	

	void getWilsonAngles () {
		double wilson_angle;
		double phijk; 
		
		arma::vec e1;
		arma::vec e2;
		arma::vec e3;

		std::cout << "Wilson angles [deg]" << std::endl;
                for (int center = 0; center < atom_vec.size(); center++) {
		std::vector<double> bondedto;
			for (int bond = 0; bond < bond_vec.size(); bond++) {
				if (bond_vec[bond].atom1 == center) {
					bondedto.push_back(bond_vec[bond].atom2);
				} else if (bond_vec[bond].atom2 == center) {
					bondedto.push_back(bond_vec[bond].atom1);
				}
			}
	
			if (bondedto.size() >= 3) { 

                        	for (int i = 0; i < bondedto.size(); i++) {
					for (int j = i + 1; j < bondedto.size(); j++) {
						for(int k = 0; k < bondedto.size(); k++) {  											// identify a triplet of bonds
							if (k != i && k != j) {
				 				e1 = getUnitVector(xyz_mat.row(center).t(), xyz_mat.row(bondedto[i]).t());  
								e2 = getUnitVector(xyz_mat.row(center).t(), xyz_mat.row(bondedto[j]).t());
								e3 = getUnitVector(xyz_mat.row(center).t(), xyz_mat.row(bondedto[k]).t());

								phijk = acos(arma::dot(e1, e2)); 
								wilson_angle = asin(arma::dot(std::pow(sin(phijk), -1) * arma::cross(e1, e2), e3));

								wilson_angles.push_back(wilson_angle);
								std::cout << bond_vec[center].atom2 << "-" << bond_vec[center].atom1 << "-" << bond_vec[i].atom2 << "-" << bond_vec[k].atom2  << " " << wilson_angle * (180.0/acos(-1.0)) << std::endl;
							}
						}
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
							argument = std::max(-1.0, std::min(1.0, argument)); 							// Had to clamp the argument because nan values appearing instead of 180deg due to the argument being -1.00000000001
							torsional_angle = acos(argument); 
	
		
							if (sin(phijk) >= bond_angle_tol or sin(phjkl) >= bond_angle_tol) { 					// Skip unphysical torsions	
								torsional_angles.push_back(torsional_angle);
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
                        for (int j = i+1; j < atom_vec.size(); j++) {                                                                                           // Loop over all unique pairs
                                if (!arma::any(neighbor_vec[i].second == j)) {
					for (int l : neighbor_vec[i].second) {
                                        	test2 = arma::join_vert(test2, neighbor_vec[l].second);
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

	// Construct molecule: 
	Molecule mol;
	mol.nr_atoms = input_mat.n_rows; 
	mol.atom_vec = input_mat.col(0);
	mol.enumerated_atom_vec = enumerateVector(input_mat.col(0));
	mol.xyz_mat  = input_mat.cols(1,3); 
	mol.makeNeighborList(); 
        mol.makeChemicalBonds();
	mol.getBondAngles();
	mol.getWilsonAngles();
	mol.getTorsionalAngles();
	mol.getInteractingAtoms();
return 0;
}



// Torsional debugging
/* 
std::cout << "\n--- Debugging Torsion Analysis ---" << std::endl;
std::cout << "Bond index: " << bond << std::endl;
std::cout << "Atom indices: atom1 = " << bond_vec[bond].atom1
          << ", atom2 = " << bond_vec[bond].atom2 << std::endl;
std::cout << "Atom types: atom1 = " << atom_vec[bond_vec[bond].atom1]
          << ", atom2 = " << atom_vec[bond_vec[bond].atom2] << std::endl;

std::cout << "Neighbors of atom1 (" << bond_vec[bond].atom1 << "): ";
for (int n : neighbor_vec[bond_vec[bond].atom1].second) {
    std::cout << n << " ";
}
std::cout << "\nNeighbors of atom2 (" << bond_vec[bond].atom2 << "): ";
for (int n : neighbor_vec[bond_vec[bond].atom2].second) {
    std::cout << n << " ";
}

std::cout << "\nTYPE OF Neighbors of atom1 (" << bond_vec[bond].atom1 << "): ";
for (int n : neighbor_vec[bond_vec[bond].atom1].second) {
    std::cout << atom_vec[n] << " ";
}
std::cout << std::endl;

std::cout << "\nTYPES OFNeighbors of atom2 (" << bond_vec[bond].atom2 << "): ";
for (int n : neighbor_vec[bond_vec[bond].atom2].second) {
    std::cout << atom_vec[n] << " ";
}
std::cout << std::endl;

std::cout << "Selected neighbors: neighbor1 = " << neighbor1
          << ", neighbor2 = " << neighbor2 << std::endl;
std::cout << "-----------------------------------\n" << std::endl;
*/


