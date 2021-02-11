#ifndef CSVPARSER_H
#define CSVPARSER_H

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

using Column =  pair<string,vector<string>>;

class CsvParser { // Class used to read csv files and store column information std::pairs (pair having string for column name and vector of strings for values)
public:
  CsvParser(string file_name);
  vector<string> GetColumn(int col);
private:
  string CsvName;
  vector<Column> Columns;
 
};





#endif 
