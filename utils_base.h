/**
 * @file   utils_base.h
 * @brief  Util functions for xml inputs. 
 */

#ifndef UTILS_BASE_H
#define UTILS_BASE_H

#include <string>
#include <boost/algorithm/string.hpp>
#include <exception>
#include <iostream>
#include <iterator>


  //-------------------------------------------------------------------------
  // String conversion utilities

  /**
   * @brief Convert a string to a number for some types
   * 
   * @param str The string we want to convert to a number
   * @retval val The string into a number
   */
  /*!The functions were not compiling. To fix it I put them as
   * inline functions. Actually I think that they should be inline function
   * because the specification of the function is so short and they are
   * good candidate for inline functions.
   */
  template<typename Number>
  Number string_to_num(
    std::string str);

  /**
   * @brief Convert a string to a @c unsigned @c char
   * 
   * @param str The string we want to convert to a number
   * @retval val The string into a number
   */
  template<>
  inline  unsigned char string_to_num<unsigned char>(
    std::string str)
  {
  return atoi(str.c_str());
  }

  /**
   * @brief Convert a string to a @c unsigned @c int 
   * 
   * @param str The string we want to convert to a number
   * @retval val The string into a number
   */
  template<>
  inline  unsigned int string_to_num<unsigned int>(
    std::string str)
  {
  return atoi(str.c_str());
  }

  /**
   * @brief Convert a string to a @c int 
   * 
   * @param str The string we want to convert to a number
   * @retval val The string into a number
   */
  template<>
  inline  int string_to_num<int>(
    std::string str)
  {
  return atoi(str.c_str());
  }

  /**
   * @brief Convert a string to a @c double
   * 
   * @param str The string we want to convert to a number
   * @retval val The string into a number
   */
  template<>
  inline  double string_to_num<double>(
    std::string str)
  {
  return atof(str.c_str());
  }

  /*!
   * @brief Convert a number of different types to a string
   * 
   * The functions were not compiling. To fix it I put them as
   * inline functions. Actually I think that they should be inline function
   * because the specification of the function is so short and they are
   * good candidate for inline functions.
   * 
   * @param val The number we want to convert to a string
   * @retval str The string we want to convert to a number
   */
  template<typename Number>
  std::string num_to_string(Number val)
  {
    // std::to_string is a C++11 standard
    return std::to_string(val);

  }


//----------------------------------------------------------------------


  template<typename Number>
  void
  string_to_vector(std::string & in, std::vector<Number> & out)
  {
    std::vector<std::string> strs;
    boost::trim(in);
    boost::split(strs,in,boost::is_any_of(", "), boost::token_compress_on);
    for (unsigned int i = 0; i<strs.size(); ++i)
      out.push_back(string_to_num<Number>(strs[i]));
  }
  
  template<typename Number>
  void
  string_to_vector(std::string & in, std::vector<std::vector<Number>> & out)
  {
    boost::trim_if(in, boost::is_any_of("\n "));
    boost::trim_right_if(in, boost::is_any_of(";"));
    std::vector<std::string> strs;
    boost::split(strs,in,boost::is_any_of(";"), boost::token_compress_on);
    out.resize(strs.size());
    for (unsigned int i = 0; i < strs.size(); ++i)
    {
      string_to_vector(strs[i], out[i]);
    }
  }

  template<typename Number>
  void
  vector_to_string(std::vector<Number> & in, std::string & out)
  {
    // We should check that the vector is not empty
    assert(in.size() != 0);

    // We put the numbers inside the string    
    std::stringstream ss;
    std::copy(in.begin(), in.end(), std::ostream_iterator<Number>(ss, " "));
    out = ss.str();
    out = out.substr(0, out.length()-1);
  }

  // We only use the indent for the pretty output
  template<typename Number>
  void
  vector_to_string(std::vector<std::vector<Number>> & in,
                   std::string & out,
                   const std::string & indent = std::string(""),
                   const unsigned int n_indent = 0)
  {
    out.clear(); out += "\n";
    std::string bin;
    for (unsigned int i = 0; i < in.size(); ++i)
    {
      vector_to_string(in[i], bin);
      for (unsigned int j = 0; j < n_indent; ++j)
        out += indent;
      out += bin + ";\n"; // we finish with a 
    }
    for (unsigned int j = 0; j < n_indent-1; ++j)
      out += indent;
  }

//----------------------------------------------------------------------

#endif

