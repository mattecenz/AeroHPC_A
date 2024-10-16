#ifndef AEROHPC_A_VTK_H
#define AEROHPC_A_VTK_H

#include <vector>
#include <array>
#include <algorithm>

// Preprocessing macros useful to checks whether the data type is supported ---- //
template<typename T>
constexpr std::string getTypeName() {
    if constexpr (std::is_same<T, double>::value) return "double";
    else if constexpr (std::is_same<T, float>::value) return "float";
    else if constexpr (std::is_same<T, int>::value) return "int";
    else static_assert(false, "Unknown data type");
}
// ----------------------------------------------------------------------------- //

/**
  * Util function that allows to swap endian definition of binary values representation
  */
template<typename T>
void swap_endian(T &var) {
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

/**
 * Stores .vtk generic data
 */
class DataSection {
public:


    /**
     * Construct a named data field with uninitialized data
     */
    DataSection(std::string section_name) : _section_name(std::move(section_name)) {
        // ensures name is underscore spaced
        std::replace(_section_name.begin(),
                     _section_name.end(),
                     ' ', '_');
    }

    /**
     * Prints .vtk binary formatted content of the class
     */
    virtual std::ostream &operator>>(std::ostream &output) = 0;

    /**
     * Returns the number of tuples contained into the data section
     */
    virtual size_t getSize() = 0;

protected:

    /**
     * Name of the data field
     */
    std::string _section_name;
};

/**
 * Stores and export .vtk formatted file
 */
class VTKFile {
public:
    /**
     * Build a class with the given file description
     */
    VTKFile(std::string description = "") : _description(description) {
        // ensure description maximum lenght
        if (_description.length() > 255)
            _description = _description.substr(0, 255) + "\n";
    }

    /**
     * Construct a class with the given file description
     * Initialize the dataset with the given physical dimension and node number
     */
    VTKFile(std::array<double, 3> physical_dimension, std::array<size_t, 3> nodes, std::string description = "")
        : VTKFile(description) {
        setupDataset(physical_dimension, nodes);
    }

    /**
     * Initialize the dataset with the given physical dimension and node number
     */
    void setupDataset(std::array<double, 3> physical_dimension, std::array<size_t, 3> nodes) {
        _nodes = nodes;
        _spacing = {physical_dimension[0] / nodes[0],
                    physical_dimension[1] / nodes[1],
                    physical_dimension[2] / nodes[2]};
        _data_size = nodes[0] * nodes[1] * nodes[2];

        _valid_dataset = true;
    }

    /**
     * Append a data section to the file if data size match the dataset size
     */
    VTKFile &operator<<(DataSection &dataSection) {
        if (_valid_dataset) {
            if (dataSection.getSize() == _data_size) {
                _data_sections.push_back(&dataSection);
                _valid_data = true;
            } else {
                fprintf(stderr, "VTKFile: Tried to insert data that doesn't match dataset");
            }
        } else {
            fprintf(stderr, "VTKFile: Tried to insert data before setup dataset");
        }
        return *this;
    }

    /**
     * Print the content of the class into a file with the given name
     */
    bool printFile(std::string file_name) {
        // ensures name is underscore spaced
        std::replace(file_name.begin(),
                     file_name.end(),
                     ' ', '_');

        // ensure file extension
        if (!file_name.ends_with(".vtk"))
            file_name += ".vtk";


        if (!_valid_dataset || !_valid_data) {
            fprintf(stderr, "VTKFile: {%s} Tried to print file before setup", file_name.c_str());
        }

        std::ofstream file;
        file.open(file_name, std::ios::out | std::ios::binary);

        printHeader(file);

        printDataset(file);

        file << "POINT_DATA " << _data_size << "\n";
        for (DataSection *data: _data_sections) {
            *data >> file;
        }

        return true;
    }

private:
    /**
     * Description of the file
     */
    std::string _description;
    /**
     * Nodes number of dataset
     */
    std::array<size_t, 3> _nodes;
    /**
     * Node spacing of the dataset
     */
    std::array<double, 3> _spacing;
    /**
     * Size of dataset (total number of nodes)
     */
    size_t _data_size;
    /**
     * Data sections appended to this class
     */
    std::vector<DataSection *> _data_sections;
    /**
     * If dataset has been correctly setup
     */
    bool _valid_dataset = false;
    /**
     * If data has been correctly setup
     */
    bool _valid_data = false;

    /**
     * Prints the header representation into the given output stream
     */
    std::ostream &printHeader(std::ostream &output) {
        output << "# vtk DataFile Version 2.0\n";
        output << _description << "\n";
        output << "BINARY\n";
        return output;
    }

    /**
     * Prints the dataset representation into the given output stream
     */
    std::ostream &printDataset(std::ostream &output) {
        output << "DATASET STRUCTURED_POINTS\n";
        output << "DIMENSIONS " << _nodes[0] << " " << _nodes[1] << " " << _nodes[2] << "\n";
        output << "SPACING " << _spacing[0] << " " << _spacing[1] << " " << _spacing[2] << "\n";
        output << "ORIGIN 0.0 0.0 0.0\n";
        return output;
    }
};


/**
 * Stores .vtk scalars array of values
 */
template<typename Value_Type, size_t Value_Size>
class ScalarsSection : public DataSection {

public:
    /**
     * Data entry value type
     */
    typedef typename std::array<Value_Type, Value_Size> entry_T;

    /**
      * Construct a named data field with uninitialized data
      */
    ScalarsSection(std::string section_name) : DataSection(section_name) {}

    /**
     * Construct a named data field with initialized data
     */
    ScalarsSection(std::string section_name, std::vector<entry_T> &values)
            : DataSection(section_name), _values(values) {}

    // Overrided functions -------------------------------------------------------- //
    std::ostream &operator>>(std::ostream &output) override {
        // standard data field header
        output << "SCALARS " << _section_name << " " << _type_name << " " << Value_Size << "\n";
        output << "LOOKUP_TABLE default\n";

        // binary of all values
        for (auto val: this->_values) {
            for (int i = 0; i < Value_Size; ++i) {
                // assign temp value and swap endian representation
                Value_Type t = val[i];
                swap_endian(t);
                // print binary
                output.write(reinterpret_cast<char *>(&t), sizeof(Value_Type));
            }
        }
        output << "\n";

        return output;
    }

    size_t getSize() override { return _values.size(); }
    // ---------------------------------------------------------------------------- //
    
    /**
     * Setup the section values
     */
    ScalarsSection &operator<<(std::vector<entry_T> &values) {
        _values = values;
        return *this;
    }

protected:

    /**
     * String representation of the data-field-entry-values type
     */
    const std::string _type_name = getTypeName<Value_Type>();

    /**
     * Data entry list
     */
    std::vector<entry_T> _values;
};


#endif //AEROHPC_A_VTK_H
