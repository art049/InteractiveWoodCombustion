// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
// 
// Permission granted to use this code only for teaching projects and 
// private practice.
//
// Do not distribute this code outside the teaching assignements.                                                                           
// All rights reserved.                                                       
// --------------------------------------------------------------------------

#pragma once

#include <string>

class Exception {
public:
  inline Exception (const std::string & msg) : m_msg ("Exception: " + msg) {}
  inline const std::string & msg () const { return m_msg; }
protected:
  std::string m_msg;
};