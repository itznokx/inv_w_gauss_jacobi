#include <cmath>
#include <iostream>
#include <vector>


std::vector<double> gauss_jacobi(std::vector<std::vector<double>> A, std::vector<double> b, double epsilon, size_t max_iter);


void print_vector(std::vector<double> input) {
  for (int i = 0; i < input.size(); i++) {
    std::cout << input[i] << ",";
  }
  std::cout << "\n";
}

double calculate_norm(std::vector<double> previous_guess, std::vector<double> current_guess) {
  double norm_num = 0;
  double norm_dem = 0;

  for (int i = 0; i < previous_guess.size(); i++) {
    double temp = std::abs(current_guess[i] - previous_guess[i]);
    
    if (temp > norm_num) {
      norm_num = temp;
    }

    if (std::abs(current_guess[i]) > norm_dem) {
      norm_dem = std::abs(current_guess[i]);
    }
  }

  return norm_num / norm_dem;
}

std::vector<double> gauss_jacobi(std::vector<std::vector<double>> A, std::vector<double> b, double epsilon, size_t max_iter) {
  std::vector<double> previous_guess;
  std::vector<double> current_guess;
  int length = A[0].size();
  
  // Populate vectors
  for (int i = 0; i < length; i++) {
    previous_guess.push_back(0);
    current_guess.push_back(0);
  }
  
  int step = 0;
  do {
    for (int i = 0; i < length; i++) {
      double ratio = 1 / A[i][i];
      
      previous_guess[i] = current_guess[i];
      current_guess[i] = 0;

      for (int j = 0; j < length; j++) {
        if (i != j) {
          current_guess[i] += (-1) * A[i][j] * previous_guess[j];
        }
      }

      current_guess[i] = (current_guess[i] + b[i]) * ratio;
    }
    step++;
  }
  while (step != max_iter && calculate_norm(previous_guess, current_guess) > epsilon);

  return current_guess;
}
