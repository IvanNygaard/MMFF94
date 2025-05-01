// Function used to compute the distance between two atoms: 

double euclideanDistance(arma::vec vec1, arma::vec vec2) {
	double distance; 
	distance = std::pow(std::pow(vec1(0) - vec2(0), 2) + std::pow(vec1(1) - vec2(1), 2) + std::pow(vec1(2) - vec2(2), 2), 0.5);
	
	return distance; 
}


arma::vec getUnitVector(arma::vec vec1, arma::vec vec2) { 
	arma::vec unitvector;
	double x_comp;
	double y_comp;
	double z_comp; 
	double distance;

	distance = euclideanDistance(vec1, vec2); 

	x_comp = -(vec1(0) - vec2(0)) * std::pow(distance, -1); 
	y_comp = -(vec1(1) - vec2(1)) * std::pow(distance, -1); 
	z_comp = -(vec1(2) - vec2(2)) * std::pow(distance, -1); 

	unitvector = {x_comp, y_comp, z_comp}; 

	return unitvector; 
}
