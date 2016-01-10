#include "parse_conf.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <unistd.h>

#define MAX_LEN_NAME 512

void help_msg(char** argv)
{
	std::cout << "Use: " << argv[0] << " [OPTIONS]" << std::endl;
	std::cout << "OPTIONS:" << std::endl;
	std::cout << "\t-h\t\t\tFor read this text" << std::endl;
	std::cout << "\t-conf <file_name>\tFor reading parameters from <file_name> (default: ./input.conf)\n" << std::endl;
}

std::map<std::string, std::string> parse_file(char * file_name)
{
	std::map<std::string, std::string> conf_map;
	std::ifstream file (file_name);
	if (!file.is_open())
	{
		std::cerr << "Error! Cannot open conf file!" << std::endl;
		exit(103);
	}
	
	char * line = new char [MAX_LEN_NAME];
	while (!file.eof())
	{
		file.getline(line, MAX_LEN_NAME);
		std:: string line_str;
		for (int i(0); i < strlen(line); ++i)
		{
			if (line[i] != ' ' && line[i] != '\0' && line[i] != '\n' && line[i] != '\t')
			{
				line_str += line[i];
			}
		}
		int find_index = line_str.find("#");
		if (find_index > -1)
		{
			line_str.erase(line_str.begin() + find_index, line_str.end());
		}
		if (line_str.length() > 0)
		{
			find_index = line_str.find("=");
			if (find_index < 0)
			{
				std::cerr << "Error! Incorrect format of conf file!" << std::endl;
				exit(104);
			}
			std::string key = line_str.substr(0, find_index);
			std::string strval = line_str.substr(find_index+1, line_str.length());
			conf_map[key] = strval;
		}
	}
	
	delete line;
	file.close();
	return conf_map;
}

std::map<std::string, std::string> parse_params(int argc, char** argv)
{
	std::map<std::string, std::string> conf_map;
	char * file_name = new char [MAX_LEN_NAME];
	char * app_dir = new char [MAX_LEN_NAME];
	int slash_index = 0;
	for (int i (strlen(argv[0]) - 1); i > 0; i--)
	{
		if (argv[0][i] == '/' || argv[0][i] == '\\')
		{
			slash_index = i+1;
			break;
		}
	}
	for (int i(0); i < slash_index; ++i)
	{
		app_dir[i] = argv[0][i];
		app_dir[i+1] = '\0';
	}
	std::string input_path = app_dir;
	if (argc == 1)
	{
		strcat(app_dir, "input.conf");
		strcpy(file_name, app_dir);
	}
	else
	{
		if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-help"))
		{
			help_msg(argv);
			exit(0);
		}
		else if (!strcmp(argv[1], "-conf"))
		{
			if (argc < 3)
			{
				std::cerr << "Error! You didn't pass file name!" << std::endl;
				exit(101);
			}
			else
			{
				slash_index = 0;
				char * in_dir = new char [MAX_LEN_NAME];
				for (int i (strlen(argv[2]) - 1); i > 0; i--)
				{
					if (argv[2][i] == '/' || argv[2][i] == '\\')
					{
						slash_index = i+1;
						break;
					}
				}
				for (int i(0); i < slash_index; ++i)
				{
					in_dir[i] = argv[2][i];
					in_dir[i+1] = '\0';
				}
				input_path = in_dir;
				delete in_dir;
				strcpy(file_name, argv[2]);
			}
		}
		else
		{
			std::cerr << "Error! Unknown argument! Use: \"" << argv[0] << " -h\" for help." << std::endl;
			exit(102);
		}
	}
	
	char * buf = new char [MAX_LEN_NAME];
    int count = readlink(file_name, buf, sizeof(char)*MAX_LEN_NAME);
    if (count >= 0) {
        buf[count] = '\0';
        strcpy(file_name, buf);
    }
	delete buf;
	
	conf_map = parse_file(file_name);
	conf_map["input_path"] = input_path;
	delete file_name;
	delete app_dir;
	
	return conf_map;
}

double ** parse_file_matrix(char* path, char* file_name)
{
	char * fname = new char [MAX_LEN_NAME];
	if (file_name[0] == '"' || file_name[strlen(file_name)-1] == '"')
	{
		int j (0);
		for (int i(0); i < strlen(file_name); ++i)
		{
			if ((i == 0 || i == strlen(file_name)-1) && file_name[i] == '"')
				continue;
			file_name[j++] = file_name[i];
		}
		file_name[j] = '\0';
	}
	strcpy(fname, path);
	strcat(fname, file_name);
	char * buf = new char [MAX_LEN_NAME];
    int count = readlink(fname, buf, sizeof(char)*MAX_LEN_NAME);
    if (count >= 0) {
        buf[count] = '\0';
        strcpy(fname, buf);
    }
	delete buf;
	std::ifstream file (fname);
	if (!file.is_open())
	{
		std::cerr << "Error! Cannot open " << fname << "!" << std::endl;
		exit(105);
	}
	int line_count = 0;
	int num_in_line = 0;
	std::string line;
	while (!file.eof())
	{
		getline(file, line);
		if (line.length() == 0)
			continue;
		if (line_count == 0)
		{
			int found = line.find(";"); 
			while(found > -1)
			{
				num_in_line++;
				found = line.find(";", found + 1);
			}
			num_in_line++;
		}
		line_count++;
	}
	file.clear();
	file.seekg (0, file.beg);
	if (line_count == 0 || num_in_line == 0)
	{
		std::cerr << "Error! Incorrect data in file " << file_name << "!" << std::endl;
		exit(106);
	}
	double ** matrix = new double * [line_count];
	for (int i(0); i < line_count; ++i)
	{
		matrix[i] = new double [num_in_line];
		getline(file, line);
		int count_nums = 0;
		int old_found = 0;
		int found;
		do
		{
			found = line.find(";", old_found + 1);
			if (found < 0)
			{
				std::cerr << "Error! Incorrect matrix size in file " << file_name << "!" << std::endl;
				exit(107);
			}
			std::string number = line.substr(old_found + 1, found);
			old_found = found;
			matrix[i][count_nums++] = atof(number.c_str());
		}
		while (count_nums < num_in_line - 1);
		std::string number = line.substr(old_found + 1, line.length());
		matrix[i][count_nums++] = atof(number.c_str());
	}
	file.close();
	return matrix;
}
