//
//  ArgumentParser.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/25/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ArgumentParser_hpp
#define ArgumentParser_hpp


#include <iostream>
#include <cstring>
#include <map>
#include <sstream>
#include <ctype.h>

// ==== helper methods : convert string to type (and specializations for string and bool) ===== //
template<typename T>
inline T convertStrToType( const std::string & s )
{
    std::stringstream convert(s);
    
    T value;
    if(convert >> value)
        return value;
    else
    {
        std::cout << "argument value for type is invalid : \t " << s << std::endl;
        std::exit(1);
    }
    return value;
}

template<>
inline bool convertStrToType( const std::string & s )
{
    std::stringstream convert(s);
    
    bool value;
    if(convert >> std::boolalpha >> value)
        return value;
    else
    {
        std::cout << "argument value for bool type is invalid : \t " << s << std::endl;
        std::cout << "use true/false only" << std::endl;
        std::exit(1);
    }
    return value;
}

template<>
inline std::string convertStrToType(const std::string & s)
{
    return s;
}

// ==== helper methods : convert type to string (and specializations for string and bool) ===== //
template<typename T>
inline std::string convertTypeToStr(const T & t)
{
    return std::to_string(t);
}

template<>
inline std::string convertTypeToStr(const bool & t)
{
    return (t ? "true" : "false");
}

template<>
inline std::string convertTypeToStr(const std::string & t)
{
    return t;
}

// ==== base class Parser ==== //
class Parser
{
protected:
    std::map<std::string,std::string> arguments_map;
    bool verbose, saveDefaults;
    
    template<typename T>
    T parse_(const std::string argkey, const bool strict, const T defval = T()) const
    {
        std::map<std::string,std::string>::const_iterator it = arguments_map.find(argkey);
        
        if (it == arguments_map.end())
        {
            if (verbose) printf("%s is empty\n", argkey.data());
            if(strict)
            {
                helpers::catastrophe("Mandatory command line option is not given. Argument name: "+argkey, __FILE__, __LINE__);
            }
            if(saveDefaults)
            {
                FILE * f = getFileHandle();
                if(f!=nullptr)
                {
                    fprintf(f, "%s %s ", argkey.data(), convertTypeToStr(defval).data());
                    fclose(f);
                }
            }
            return defval;
        }
        else
        {
            if (verbose) printf("Found the value for key %s as %s\n", argkey.data(), it->second.data());
            return convertStrToType<T>(it->second);
        }
    }
    
    FILE * getFileHandle() const
    {
        const std::string filepath = "argumentparser.log";
        FILE * f = fopen(filepath.c_str(), "a");
        if (f == nullptr)
        {
            printf("Can not open file %s.\n", filepath.data());
            return nullptr;
        }
        return f;
    }
    
public:
    template<typename T>
    T parse(const std::string arg) const
    {
        return parse_<T>(arg, true);
    }
    
    template<typename T>
    T parse(const std::string arg, const T defval) const
    {
        return parse_<T>(arg, false, defval);
    }
    
    Parser(): arguments_map(), verbose(false), saveDefaults(false)
    {
    }
    
    void set_verbosity(const bool verbosity)
    {
        verbose = verbosity;
    }
    
    void save_defaults()
    {
        saveDefaults = true;
    }
    
    void print_options() const
    {
        for(const auto & it : arguments_map)
            printf("%s %s ",it.first.data(), it.second.data());
        printf("\n");
    }
    
    void save_options()
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        for(const auto & it : arguments_map)
            fprintf(f, "%s %s ", it.first.data(), it.second.data());
        fclose(f);
    }
    
    void finalize() const
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        fprintf(f, "\n");
        fclose(f);
    }
};

class ArgumentParser : public Parser
{
private:
    
    const int parse_argc;
    const char** parse_argv;
    
public:
    
    ArgumentParser(const int argc, const char ** argv) : Parser(),parse_argc(argc), parse_argv(argv)
    {
        // loop over the arguments
        for (int i=1; i<argc; i++)
        {
            // check all exceptions: we are very strict! no value-less arguments can be given
            if(argv[i][0] != '-')
                helpers::catastrophe("found an unexpected command-line entry : "+std::string(argv[i]), __FILE__, __LINE__);
            
            if(i == argc-1)
                helpers::catastrophe("the last command-line argument is without value : "+std::string(argv[i]), __FILE__, __LINE__);
            
            if( (argv[i+1][0] == '-') and not std::isdigit(argv[i+1][1]) )
                helpers::catastrophe("found a command-line argument without value : "+std::string(argv[i]), __FILE__, __LINE__);
            
            if(arguments_map.find(argv[i]) != arguments_map.end())
                helpers::catastrophe("found a multiply-defined command-line argument : "+std::string(argv[i]), __FILE__, __LINE__);
            
            arguments_map[argv[i]] = argv[i+1];
            i += 1;
        }
        
        if(verbose) printf("found %ld arguments of %d\n",arguments_map.size(),argc);
    }
    
    int getargc() const
    {
        return parse_argc;
    }
    
    const char** getargv() const
    {
        return parse_argv;
    }
};


class StringParser: public Parser
{
    // Can call this parser with a string of arguments like
    // std::string testString = "arg1=val1 arg2=val2 arg3=val3";
    // StringParser parser(testString);
    // std::cout << "arg1 = " << parser("-arg1").asString() << std::endl;
    
protected:
    
    // trim methods from stackoverflow
    
    // trim from start
    inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }
    
    // trim from end
    inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }
    
    // trim from both ends
    inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
    }
    
public:
    
    StringParser(std::string & input):
    Parser()
    {
        std::string key,value;
        std::istringstream iss(input);
        
        while( std::getline(iss, key, '=') )
        {
            if( std::getline(iss, value, ' ') )
            {
                // add "-" because then we can use the same code for parsing factory as command lines
                arguments_map["-"+trim(key)] = trim(value);
            }
        }
        
        set_verbosity(false);
    }
};

#endif /* ArgumentParser_hpp */
