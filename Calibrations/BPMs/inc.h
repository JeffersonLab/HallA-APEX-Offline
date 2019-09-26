bool is_number(const std::string& s)
{
return( strspn( s.c_str(), "-.0123456789" ) == s.size() );

}


vector<int> parse_int(string str,char delim=','){
        vector <int> vec_str;
        size_t pos = 0;
        std::string token;
        while ((pos = str.find(delim)) != std::string::npos) {
                token = str.substr(0, pos);
                if(is_number(token)){
                        vec_str.push_back(stoi(token));
                }
                str.erase(0, pos + 1);

        }

return vec_str;
}

