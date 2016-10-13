/*
 * Arguments 2.0 Beta - A Command Line Processing Library
 * Copyright (C) 2000, 2001 Jared Davis
 *
 * This program is free software; you can redistribute it and/or modif2 it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "Arguments.h"
#include <algorithm>

std::vector<std::string> Arguments::s_spaceArgs;

//
// ExplodeString function
// ----------------------
//   Explodes a string into tokens, which it adds to a vector.
//   (This function is taken from the AFA library, also by Jared Davis)
//
void Arguments::ExplodeString(const std::string& str, std::vector<std::string>& tokens, char delimiter)
throw(std::bad_alloc) {
    std::string::size_type next, prev = 0;

    do {
        next = str.find(delimiter, prev);

        tokens.push_back( str.substr(prev, next - prev) );

        prev = next + 1;

    } while (next != std::string::npos);
}

//
// findArgument
// ------------
//   Performs a binary search on the argument vector and returns the integer
//   corresponding to the argument, if it exists.  Otherwise returns -1.
//   arg is a pipe-delimited vector of acceptable arguments.
//
int Arguments::findArgument(const std::string& arg) const throw() {
    if (d_arguments.size() == 0) {
        return -1;
    }

    std::vector<std::string> tokens;
    ExplodeString(arg, tokens, '|');

    for (unsigned int i = 0; i < tokens.size(); ++i) {
        int low = 0;
        int high = d_arguments.size() - 1;

        while (low <= high) {
            unsigned int mid = (low + high) / 2;

            if (tokens[i] < d_arguments[mid].first) {
                high = mid - 1;
            }

            else if (tokens[i] > d_arguments[mid].first) {
                low = mid + 1;
            }

            else {
                return mid;
            }
        }
    }

    return -1;
}


//
// Arguments::setArgumentsWithSpaces
// ---------------------------------
//   Adds all tokens of the pipe-delimited string, args, to the
//   spaces vector.  keeps the vector sorted.
//
#include <iostream>
void Arguments::setArgumentsWithSpaces(const std::string& args)
throw(std::bad_alloc) {
    ExplodeString(args, s_spaceArgs, '|');
    std::sort(s_spaceArgs.begin(), s_spaceArgs.end());
}


Arguments::Arguments(int argc, const char** argv) throw(std::bad_alloc) {
    // start at 1 to skip program's name.
    for (int i = 1; i < argc; ++i) {
        // Extract this token.
        std::string this_arg = argv[i];

        // Try to find a colon or an equals sign.
        std::string::size_type colon = this_arg.find(':');
        if (colon == std::string::npos) {
            colon = this_arg.find('=');
        }

        // Add arg/val pair if we have a colon or equals sign.
        if (colon != std::string::npos) {
            d_arguments.push_back(std::make_pair(
                                      this_arg.substr(0, colon),
                                      std::make_pair(this_arg.substr(colon + 1), true)
                                  ));

            continue;
        }

        // If we get here, we didn't find a colon or an equals, we need to
        // consult the spaces vector.
        if (i < (argc - 1)) {
            if (std::binary_search(s_spaceArgs.begin(), s_spaceArgs.end(), this_arg)) {
                // if parameter value starts with "-" consider it as an option
                if (argv[i+1][0]=='-') {
                    d_arguments.push_back(std::make_pair(
                                              this_arg,
                                              std::make_pair(std::string(""), false)));
                    continue;   // i is not incrented !!
                } else {
                    d_arguments.push_back(std::make_pair(
                                              this_arg,
                                              std::make_pair(argv[++i], true)
                                          ));
                    continue;
                }
            }
        }

        // If we get here, we either don't have a space, or we are at the last
        // argument that we were passed, so we can't possibly have a value for
        // this argument. Throw it into the map with an empty value.
        d_arguments.push_back(std::make_pair(
                                  this_arg,
                                  std::make_pair(std::string(""), false)
                              ));
    }

    // Now let's sort the argument vector so it's fast to find them.
    std::sort(d_arguments.begin(), d_arguments.end());
}


//
// Arguments::size
// ---------------
//   Returns the number of command-line arguments that your program was passed.
//   This is not the same as argc: it does not include your program's name, nor
//   does it include the values passed to your command-line arguments.
//
unsigned int Arguments::size() const throw() {
    return d_arguments.size();
}

//
// Arguments::has
// --------------
//   Cycle through the tokens of argument and return true if any of them
//   are found.
//
bool Arguments::has(const std::string& argument) const throw(std::bad_alloc) {
    return findArgument(argument) != -1;
}




/*
#ifndef NDEBUG

	#include <iostream>

	void Arguments::debug() const
	{
		std::cout << "Arguments Debug Information-----\n" << std::endl;
		std::cout << "Number of space-taking arguments: " << s_spaceArgs.size() << "\n\t";
		copy(s_spaceArgs.begin(), s_spaceArgs.end(), ostreaiterator_<string>(std::cout, "\n\t"));

		std::cout << "Arguments Map: " << size() << " Elements. " << std::endl;
		for(unsigned int i = 0;i < size();++i)
		{
			std::cout << "\t" << d_arguments[i].first << " = " << d_arguments[i].second.first << "   ("
				<< d_arguments[i].second.second << ")\n";
		}
		std::cout << "\n-----End of Arguments Debug Information" << std::endl;
	}

#endif
*/
