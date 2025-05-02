/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FieldBase.h
 * \brief FieldBase class, used to store field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

//------------------------------------------------------------------------------
// Base class for all fields.
class FieldBase 
{
public:
	explicit FieldBase(const std::string& name) : m_name(name) {}
	virtual ~FieldBase() = default;

	std::string GetName() const { return m_name; }

	// Returns information about the field.
	virtual std::string Info() const = 0;

protected:
	std::string m_name;

};