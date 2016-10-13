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

#ifndef INCLUDED_OPENDELI_ARGUMENTS_H
#define INCLUDED_OPENDELI_ARGUMENTS_H

#if (defined(_WIN32) && defined(_DEBUG))
// In visual studio, compiling this library gives you a billion and a half
// error messages that "debug information has been truncated to 255 chars"
// because the visual studio debugger can't handle symbols longer than 255
// characters.  This disables this warning so we comile cleanly.
#pragma warning( disable: 4786 )
#endif

#include <vector>
#include <string>
#include <sstream>



class Arguments {
private:

    // Vector of the arguments that take spaces.
    static std::vector<std::string> s_spaceArgs;

    // Vector of argument-value pairs.
    std::vector< std::pair<std::string, std::pair<std::string, bool> > > d_arguments;


    // explodes a delimited string into its tokens.
    static void ExplodeString(const std::string& str, std::vector<std::string>& tokens, char delimiter)
    throw(std::bad_alloc);

    // finds the index of the argument arg.
    int findArgument(const std::string& arg) const throw();

    //
    // convert(source, target)
    //   Converts a variable of type A into a variable of type B by passing
    //   it through a strstream.
    //
    template<class A, class B>
    static bool convert(const A& source, B& target) throw() {
        std::stringstream ss;
        ss.clear();
        ss << source;
        ss >> target;
        return (!ss.bad() && !ss.fail());
    }


public:

    static void setArgumentsWithSpaces(const std::string& args)
    throw(std::bad_alloc);

    Arguments(int argc, const char** argv)
    throw(std::bad_alloc);

    unsigned int size() const throw();

    bool has(const std::string& arg) const throw(std::bad_alloc);

    template<class A, class B> bool get
    (const std::string& arg, A& value, const B& default_value) const
    throw(std::bad_alloc) {
        int index = findArgument(arg);

        if (index == -1) {
            value = default_value;
            return false;
        }

        else {
            bool ret = convert(d_arguments[index].second.first, value);
            if (!ret)
                value = default_value;
            return ret;
        }
    }


#ifndef NDEBUG
    void debug() const;
#endif

};

#endif // INCLUDED_OPENDELI_ARGUMENTS_H
