#include "AlignmentAlgorithm.h"
#include <iostream>
#include <fstream>

namespace seq {

double** AlignmentAlgorithm::IUB()
{
  static double rowA[] = { 5,-4,-4,-4,1,1,1,-4,-4,-4,-1,-1,-1,-4,-2 };
  static double rowC[] = { -4,5,-4,-4,1,-4,-4,1,1,-4,-1,-1,-4,-1,-2 };
  static double rowG[] = { -4,-4,5,-4,-4,1,-4,1,-4,1,-1,-4,-1,-1,-2 };
  static double rowT[] = { -4,-4,-4,5,-4,-4,1,-4,1,1,-4,-1,-1,-1,-2 };
  static double rowM[] = { 1,1,-4,-4,-1,-2,-2,-2,-2,-4,-1,-1,-3,-3,-1 };
  static double rowR[] = { 1,-4,1,-4,-2,-1,-2,-2,-4,-2,-1,-3,-1,-3,-1 };
  static double rowW[] = { 1,-4,-4,1,-2,-2,-1,-4,-2,-2,-3,-1,-1,-3,-1 };
  static double rowS[] = { -4,1,1,-4,-2,-2,-4,-1,-2,-2,-1,-3,-3,-1,-1 };
  static double rowY[] = { -4,1,-4,1,-2,-4,-2,-2,-1,-2,-3,-1,-3,-1,-1 };
  static double rowK[] = { -4,-4,1,1,-4,-2,-2,-2,-2,-1,-3,-3,-1,-1,-1 };
  static double rowV[] = { -1,-1,-1,-4,-1,-1,-3,-1,-3,-3,-1,-2,-2,-2,-1 };
  static double rowH[] = { -1,-1,-4,-1,-1,-3,-1,-3,-1,-3,-2,-1,-2,-2,-1 };
  static double rowD[] = { -1,-4,-1,-1,-3,-1,-1,-3,-3,-1,-2,-2,-1,-2,-1 };
  static double rowB[] = { -4,-1,-1,-1,-3,-3,-3,-1,-1,-1,-2,-2,-2,-1,-1 };
  static double rowN[] = { -2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

  static double *iub[] = { rowA, rowC, rowG, rowT, rowM, rowR, rowW, rowS,
			   rowY, rowK, rowV, rowH, rowD, rowB, rowN };

  return iub;
}

double** AlignmentAlgorithm::BLOSUM30()
{
  static double rowA[] = { 4,-3,0,0,-2,0,-2,0,0,-1,1,0,-1,1,-1,1,1,1,-5,-4,-7,0,0,0,0,0 };
  static double rowC[] = { -3,17,-3,1,-3,-4,-5,-2,-3,0,-2,-1,-3,-2,-2,-2,-2,-2,-2,-6,-7,0,0,0,-2,-2 };
  static double rowD[] = { 0,-3,9,1,-5,-1,-2,-4,0,-1,-3,1,-1,-1,-1,0,-1,-2,-4,-1,-7,0,0,0,5,-1 };
  static double rowE[] = { 0,1,1,6,-4,-2,0,-3,2,-1,-1,-1,1,2,-1,0,-2,-3,-1,-2,-7,0,5,0,0,-1 };
  static double rowF[] = { -2,-3,-5,-4,10,-3,-3,0,-1,2,-2,-1,-4,-3,-1,-1,-2,1,1,3,-7,0,-4,0,-3,-1 };
  static double rowG[] = { 0,-4,-1,-2,-3,8,-3,-1,-1,-2,-2,0,-1,-2,-2,0,-2,-3,1,-3,-7,0,-2,0,0,-1 };
  static double rowH[] = { -2,-5,-2,0,-3,-3,14,-2,-2,-1,2,-1,1,0,-1,-1,-2,-3,-5,0,-7,0,0,0,-2,-1 };
  static double rowI[] = { 0,-2,-4,-3,0,-1,-2,6,-2,2,1,0,-3,-2,-3,-1,0,4,-3,-1,-7,0,-3,0,-2,0 };
  static double rowK[] = { 0,-3,0,2,-1,-1,-2,-2,4,-2,2,0,1,0,1,0,-1,-2,-2,-1,-7,0,1,0,0,0 };
  static double rowL[] = { -1,0,-1,-1,2,-2,-1,2,-2,4,2,-2,-3,-2,-2,-2,0,1,-2,3,-7,0,-1,0,-1,0 };
  static double rowM[] = { 1,-2,-3,-1,-2,-2,2,1,2,2,6,0,-4,-1,0,-2,0,0,-3,-1,-7,0,-1,0,-2,0 };
  static double rowN[] = { 0,-1,1,-1,-1,0,-1,0,0,-2,0,8,-3,-1,-2,0,1,-2,-7,-4,-7,0,-1,0,4,0 };
  static double rowP[] = { -1,-3,-1,1,-4,-1,1,-3,1,-3,-4,-3,11,0,-1,-1,0,-4,-3,-2,-7,0,0,0,-2,-1 };
  static double rowQ[] = { 1,-2,-1,2,-3,-2,0,-2,0,-2,-1,-1,0,8,3,-1,0,-3,-1,-1,-7,0,4,0,-1,0 };
  static double rowR[] = { -1,-2,-1,-1,-1,-2,-1,-3,1,-2,0,-2,-1,3,8,-1,-3,-1,0,0,-7,0,0,0,-2,-1 };
  static double rowS[] = { 1,-2,0,0,-1,0,-1,-1,0,-2,-2,0,-1,-1,-1,4,2,-1,-3,-2,-7,0,-1,0,0,0 };
  static double rowT[] = { 1,-2,-1,-2,-2,-2,-2,0,-1,0,0,1,0,0,-3,2,5,1,-5,-1,-7,0,-1,0,0,0 };
  static double rowV[] = { 1,-2,-2,-3,1,-3,-3,4,-2,1,0,-2,-4,-3,-1,-1,1,5,-3,1,-7,0,-3,0,-2,0 };
  static double rowW[] = { -5,-2,-4,-1,1,1,-5,-3,-2,-2,-3,-7,-3,-1,0,-3,-5,-3,20,5,-7,0,-1,0,-5,-2 };
  static double rowY[] = { -4,-6,-1,-2,3,-3,0,-1,-1,3,-1,-4,-2,-1,0,-2,-1,1,5,9,-7,0,-2,0,-3,-1 };
  static double rowSTP[] = { -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,1,0,-7,0,-7,-7 };
  static double rowGAP[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static double rowZ[] = { 0,0,0,5,-4,-2,0,-3,1,-1,-1,-1,0,4,0,-1,-1,-3,-1,-2,-7,0,4,0,0,0 };
  static double rowU[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static double rowB[] = { 0,-2,5,0,-3,0,-2,-2,0,-1,-2,4,-2,-1,-2,0,0,-2,-5,-3,-7,0,0,0,5,-1 };
  static double rowX[] = { 0,-2,-1,-1,-1,-1,-1,0,0,0,0,0,-1,0,-1,0,0,0,-2,-1,-7,0,0,0,-1,-1 };

  static double *mat[] = { rowA, rowC, rowD, rowE, rowF, rowG, rowH, rowI,
			   rowK, rowL, rowM, rowN, rowP, rowQ, rowR, rowS,
			   rowT, rowV, rowW, rowY, rowSTP, rowGAP,
			   rowZ, rowU, rowB, rowX };

  return mat;
}

double** AlignmentAlgorithm::parseMatrix(const std::string &filename) {
    // Open the file
    std::ifstream file(filename);

    // Check if the file was successfully opened
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return nullptr;
    }

    // Consume all lines with # at the beginning
    std::string line;
    while (std::getline(file, line)) {
        if (line[0] != '#') {
            break;
        }
    }

    // First line will tell us the order of the columns. It will also tell us the size of the array.
    std::string expectedColOrder = "ARNDCQEGHILKMFPSTWYVBJZX*_";
    int columnCounter = 0;

    for (int i = 0; i < line.size(); i++){
        if(line[i] != ' '  & line[i] != ','){
            if(line[i] != expectedColOrder[columnCounter]){
                std::cerr << "Error: Columns are not in the expected order. Received " << line[i] << " instead of " << expectedColOrder[columnCounter];
                return nullptr;
            }else{
                columnCounter++;
            }
        }
    }

    cout << "Validation passed.\n";

    double** outputTable = new double*[expectedColOrder.size()];

    int rowCounter = 0;



    // Consume the rest of the file
    while (std::getline(file, line)) {

        // Check that this row's amino acid is correctly aligned with the column order
        if(line[0] == expectedColOrder[rowCounter]){
            double* rowArray = new double[expectedColOrder.size()];
            columnCounter = 0;
            std::string value;

            // We need to account for negatives, so we consume the next character until we come across a comma
            // We can also start at position 2 since pos 0 is the AA and pos 1 is a comma
            for(int i = 2; i < line.size(); i++){
                if(line[i] != ','){
                    value += line[i];
                }else{
                    // When we reach the next comma, we can convert the stored value to a double and put it in the array
                    double valueDouble = std::stod(value.c_str());
                    rowArray[columnCounter] = valueDouble;

                    columnCounter ++;
                    value = "";
                }
            }

            outputTable[rowCounter] = rowArray;
            rowCounter++;
        }else{
            std::cerr << "Error: The rows are not in the right order. They should match the column order.\n";
            std::cerr << "Expected " << expectedColOrder[rowCounter] << " but got " << line[0] << " instead";
        }
    }


    return outputTable;
    }

};
