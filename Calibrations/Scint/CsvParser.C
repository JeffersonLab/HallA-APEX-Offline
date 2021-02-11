#include "CsvParser.h"

CsvParser::CsvParser(string file_name):
  CsvName(file_name)
  {

    ifstream infile(file_name);
    string line = "";
    
    //    vector<Column> Columns;
    
    
    int line_count = 0;
    while (getline(infile, line)){
        stringstream strstr(line);
        string word = "";
	vector<string> value_line; // values for one line
	int word_count = 0;
        while (getline(strstr,word, ',')){
	  if(line_count == 0){
	    vector<string> placeholder;
	    pair<string , vector<string>> head_col(word, placeholder);
	    Columns.push_back(head_col);
	  }
	  else{
	    Columns[word_count].second.push_back(word);
	  }
	  word_count++;
	}
	line_count++;
    }
    
    
  }


vector<string> CsvParser::GetColumn(int col){
  // get specified column
  
  vector<string> Col = Columns[col].second;
  return Col;

}
