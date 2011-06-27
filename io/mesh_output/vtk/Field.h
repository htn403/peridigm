/*
 * Field.h
 *
 *  Created on: Mar 18, 2010
 *      Author: jamitch
 */

#ifndef FIELD_H_
#define FIELD_H_
#include <stdexcept>
#include <tr1/memory>
#include <string>
#include <map>
#include "utilities/Array.h"

using std::tr1::shared_ptr;
using UTILITIES::Array;
using std::string;
using std::size_t;
using std::map;

namespace Field_ENUM {

enum Relation {
			POINT=0,
			BOND,
			RELATION_UNDEFINED
};

enum Type {
	        VOLUME=0,
			GID,
			PROC_NUM,
			WEIGHTED_VOLUME,
			DILATATION,
			DAMAGE,
			E_DP,
			E_DB,
			PLASTIC_CONSISTENCY,
			NORM_DEVIATORIC_FORCE_STATE,
			DSF,
			NUM_NEIGHBORS,
			COORDINATES,
			DISPLACEMENT,
			CURRENT_COORDINATES,
			VELOCITY,
			ACCELERATION,
			BC_MASK,
			FORCE,
			FORCE_DENSITY,
			CONTACT_FORCE,
			CONTACT_FORCE_DENSITY,
			RESIDUAL,
			BOND_DAMAGE,
			TYPE_UNDEFINED
};

enum Length {
			LENGTH_UNDEFINED=0,
			SCALAR=1,
			VECTOR2D=2,
			VECTOR3D=3
};

enum Temporal {
	SCRATCH=0,
	CONSTANT,
	TWO_STEP_INTEGRATED,
	TEMPORAL_UNDEFINED
};

enum ParallelTopology {
	        OWNED=0,
	        OVERLAP,
	        NONE,
	        PARALLEL_TOPOLOGY_UNDEFINED
};

enum Step {
			STEP_N=0,
			STEP_NP1=1,
			STEP_NONE=2,
			STEP_UNDEFINED
};

}

namespace Field_NS {

struct FieldSpec  {
private:
	Field_ENUM::Type type;
	Field_ENUM::Relation relation;
	Field_ENUM::Length length;
	Field_ENUM::ParallelTopology par_top;
	Field_ENUM::Temporal temporal;
	unsigned int id;
	string label;
	FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::ParallelTopology p, Field_ENUM::Temporal temp, const string& label);

public:
	FieldSpec();
	FieldSpec(const FieldSpec& c);
	explicit FieldSpec(Field_ENUM::Type t, Field_ENUM::Relation r,  Field_ENUM::Length len, Field_ENUM::Temporal temp, const string& label);
	const FieldSpec& operator=(const FieldSpec& c);
	bool operator == (const FieldSpec& right) const;
	bool operator != (const FieldSpec& right) const;
	bool operator < (const FieldSpec& right) const;
	std::ostream& print(std::ostream& os) const;
	const Field_ENUM::Type getType() const { return type; }
	const Field_ENUM::Relation getRelation() const { return relation; }
	const Field_ENUM::Length getLength() const { return length; }
	const Field_ENUM::ParallelTopology get_parallel_topology() const { return par_top; }
	const Field_ENUM::Temporal get_temporal() const { return temporal; }
	const unsigned int getId() const { return id; }
	const string getLabel() const { return label; }
	const FieldSpec get_overlap_spec() const;
	const FieldSpec get_override(Field_ENUM::Temporal temp) const;
};

inline std::ostream& operator<<(std::ostream& os, const FieldSpec& fs) {
	  return fs.print(os);
}
/*
 * UNDEFINED FIELDSPEC
 */
const FieldSpec FIELDSPEC_UNDEFINED(Field_ENUM::TYPE_UNDEFINED,Field_ENUM::RELATION_UNDEFINED,Field_ENUM::LENGTH_UNDEFINED,Field_ENUM::TEMPORAL_UNDEFINED,"FIELDSPEC_UNDEFINED");

template<typename T>
class Field : public FieldSpec, public Array<T> {
private:
	size_t num_points;

public:
	Field() : FieldSpec(FIELDSPEC_UNDEFINED), Array<T>(0*FIELDSPEC_UNDEFINED.getLength()), num_points(0) {}
    Field(const FieldSpec& c, size_t numPoints) : FieldSpec(c), Array<T>(numPoints*c.getLength()), num_points(numPoints) {}
    Field(const FieldSpec& c, shared_ptr<T> data, size_t numPoints) :
    	FieldSpec(c), Array<T>(numPoints*c.getLength(),data), num_points(numPoints) {}
    size_t get_num_points() const { return num_points; }
};


template<typename T>
class TemporalField {
public:

	TemporalField(const FieldSpec& c=FIELDSPEC_UNDEFINED,  size_t numPoints=0): num_points(numPoints), N(c,numPoints), NP1(c,numPoints) {}
	size_t get_num_points() const { return num_points; }
    Field<T> getField(Field_ENUM::Step step) const {
    	switch(step){
    	case Field_ENUM::STEP_N:
    		return N;
    		break;
    	case Field_ENUM::STEP_NP1:
    		return NP1;
    		break;
    	default:
    		std::string m = "TemporalField::TemporalField(Invalid value for argument Step)\n";
    		throw std::invalid_argument(m);
    	}
    }

    /*
     * This function swaps the data between N and NP1;
     * Calculations are generally performed on NP1.
     * During time integration, NP1 can be thought of
     * as a scratch array;  At some point, NP1 will satisfy
     * iteration requirements, then, this call
     * advances the time step in order to save NP1 to N.
     * Then the process begins anew.
     */
    void advanceStep() { std::swap(N,NP1); }

private:
    size_t num_points;
	Field<T> N, NP1;
};

/*
 * Constant FieldSpecs defined for use in PdITI namespace
 */

/*
* POINT SCALAR FieldSpecs
*/
const Field_NS::FieldSpec VOLUME
(Field_ENUM::VOLUME, Field_ENUM::POINT,                 Field_ENUM::SCALAR,   Field_ENUM::CONSTANT,       "Volume");

const Field_NS::FieldSpec GID
(Field_ENUM::GID,Field_ENUM::POINT,                    Field_ENUM::SCALAR,    Field_ENUM::CONSTANT,          "GID");

const Field_NS::FieldSpec PROC_NUM
( Field_ENUM::PROC_NUM, Field_ENUM::POINT,             Field_ENUM::SCALAR,    Field_ENUM::CONSTANT,     "Proc_Num");

const Field_NS::FieldSpec WEIGHTED_VOLUME
(Field_ENUM::WEIGHTED_VOLUME, Field_ENUM::POINT,       Field_ENUM::SCALAR,  Field_ENUM::CONSTANT, "Weighted_Volume");

const Field_NS::FieldSpec DILATATION
(Field_ENUM::DILATATION, Field_ENUM::POINT,            Field_ENUM::SCALAR,  Field_ENUM::TWO_STEP_INTEGRATED,     "Dilatation");

const Field_NS::FieldSpec NUM_NEIGHBORS
( Field_ENUM::NUM_NEIGHBORS, Field_ENUM::POINT,        Field_ENUM::SCALAR,  Field_ENUM::CONSTANT,  "Num_Neighbors");

const Field_NS::FieldSpec LAMBDA
(Field_ENUM::PLASTIC_CONSISTENCY, Field_ENUM::POINT,    Field_ENUM::SCALAR,   Field_ENUM::TWO_STEP_INTEGRATED,       "Lambda");

const Field_NS::FieldSpec NORM_TD
(Field_ENUM::NORM_DEVIATORIC_FORCE_STATE, Field_ENUM::POINT, Field_ENUM::SCALAR,   Field_ENUM::SCRATCH,      "Norm_td");

const Field_NS::FieldSpec DSF
(Field_ENUM::DSF,   Field_ENUM::POINT,                  Field_ENUM::SCALAR,    Field_ENUM::SCRATCH,          "DSF");

const Field_NS::FieldSpec BC_MASK
(Field_ENUM::BC_MASK,  Field_ENUM::POINT,               Field_ENUM::SCALAR,    Field_ENUM::SCRATCH,      "BC_MASK");

const Field_NS::FieldSpec DAMAGE
(Field_ENUM::DAMAGE,  Field_ENUM::POINT,               Field_ENUM::SCALAR,    Field_ENUM::SCRATCH,      "Damage");

/**
 * POINT VECTOR3D FieldSpecs
 */
const Field_NS::FieldSpec COORD3D
(Field_ENUM::COORDINATES,  Field_ENUM::POINT,           Field_ENUM::VECTOR3D,     Field_ENUM::CONSTANT,       "Coordinates");

const Field_NS::FieldSpec DISPL3D
(Field_ENUM::DISPLACEMENT,   Field_ENUM::POINT,         Field_ENUM::VECTOR3D,        Field_ENUM::TWO_STEP_INTEGRATED,   "Displacement");

const Field_NS::FieldSpec CURCOORD3D
(Field_ENUM::CURRENT_COORDINATES, Field_ENUM::POINT,    Field_ENUM::VECTOR3D,   Field_ENUM::TWO_STEP_INTEGRATED, "Current_Coordinates");

const Field_NS::FieldSpec VELOC3D
(Field_ENUM::VELOCITY, Field_ENUM::POINT,               Field_ENUM::VECTOR3D,            Field_ENUM::TWO_STEP_INTEGRATED,  "Velocity");

const Field_NS::FieldSpec ACCEL3D
(Field_ENUM::ACCELERATION,  Field_ENUM::POINT,          Field_ENUM::VECTOR3D,       Field_ENUM::TWO_STEP_INTEGRATED,   "Acceleration");

const Field_NS::FieldSpec FORCE3D
(Field_ENUM::FORCE,  Field_ENUM::POINT,                 Field_ENUM::VECTOR3D,              Field_ENUM::TWO_STEP_INTEGRATED,   "Force");

const Field_NS::FieldSpec FORCE_DENSITY3D
(Field_ENUM::FORCE_DENSITY, Field_ENUM::POINT,          Field_ENUM::VECTOR3D,       Field_ENUM::TWO_STEP_INTEGRATED,  "Force Density");

const Field_NS::FieldSpec CONTACT_FORCE3D
(Field_ENUM::CONTACT_FORCE,  Field_ENUM::POINT,         Field_ENUM::VECTOR3D,    Field_ENUM::TWO_STEP_INTEGRATED,     "Contact Force");

const Field_NS::FieldSpec CONTACT_FORCE_DENSITY3D
(Field_ENUM::CONTACT_FORCE_DENSITY, Field_ENUM::POINT,  Field_ENUM::VECTOR3D, Field_ENUM::TWO_STEP_INTEGRATED, "Contact Force Density");

const Field_NS::FieldSpec RESID3D
(Field_ENUM::RESIDUAL,  Field_ENUM::POINT,             Field_ENUM::VECTOR3D,    Field_ENUM::SCRATCH,          "Residual");

/**
 * BOND SCALAR FieldSpecs
 */
const Field_NS::FieldSpec BOND_DAMAGE
(Field_ENUM::BOND_DAMAGE, Field_ENUM::BOND,  Field_ENUM::SCALAR, Field_ENUM::SCRATCH, "Bond Damage");

const Field_NS::FieldSpec DEVIATORIC_PLASTIC_EXTENSION
(Field_ENUM::E_DP, Field_ENUM::BOND,   Field_ENUM::SCALAR, Field_ENUM::TWO_STEP_INTEGRATED, "Deviatoric_Plastic_Extension");

const Field_NS::FieldSpec DEVIATORIC_BACK_EXTENSION
(Field_ENUM::E_DB, Field_ENUM::BOND,    Field_ENUM::SCALAR,  Field_ENUM::TWO_STEP_INTEGRATED, "Deviatoric_Back_Extension");


struct FieldSpecMap {
	static std::map<string, FieldSpec> create_map() {
		std::map<string,FieldSpec> mymap;
		mymap[FIELDSPEC_UNDEFINED.getLabel()]          = FIELDSPEC_UNDEFINED;
		mymap[VOLUME.getLabel()]                       = VOLUME;
		mymap[GID.getLabel()]                          = GID;
		mymap[PROC_NUM.getLabel()]                     = PROC_NUM;
		mymap[DAMAGE.getLabel()]                       = DAMAGE;
		mymap[WEIGHTED_VOLUME.getLabel()]              = WEIGHTED_VOLUME;
		mymap[DILATATION.getLabel()]                   = DILATATION;
		mymap[NUM_NEIGHBORS.getLabel()]                = NUM_NEIGHBORS;
		mymap[LAMBDA.getLabel()]                       = LAMBDA;
		mymap[DSF.getLabel()]                          = DSF;
		mymap[COORD3D.getLabel()]                      = COORD3D;
		mymap[DISPL3D.getLabel()]                      = DISPL3D;
		mymap[CURCOORD3D.getLabel()]                   = CURCOORD3D;
		mymap[VELOC3D.getLabel()]                      = VELOC3D;
		mymap[ACCEL3D.getLabel()]                      = ACCEL3D;
		mymap[FORCE3D.getLabel()]                      = FORCE3D;
		mymap[FORCE_DENSITY3D.getLabel()]              = FORCE_DENSITY3D;
		mymap[CONTACT_FORCE3D.getLabel()]              = CONTACT_FORCE3D;
		mymap[CONTACT_FORCE_DENSITY3D.getLabel()]      = CONTACT_FORCE_DENSITY3D;
		mymap[BOND_DAMAGE.getLabel()]                  = BOND_DAMAGE;
		mymap[DEVIATORIC_PLASTIC_EXTENSION.getLabel()] = DEVIATORIC_PLASTIC_EXTENSION;
		mymap[DEVIATORIC_BACK_EXTENSION.getLabel()]    = DEVIATORIC_BACK_EXTENSION ;
		return mymap;
	};
	static const std::map<string, FieldSpec> Map;
};

} // Field_NS

#endif // FIELD_H_