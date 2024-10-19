#ifndef AEROHPC_A_VTK_H
#define AEROHPC_A_VTK_H

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <fstream>
#include <functional>
#include "Traits.hpp"

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
 * Util function that reads a value of type T from the input stream
 */
template<typename T>
static inline bool getword(std::ifstream &input, T &buff) {
    bool valid = false;
    if (input >> buff) valid = true;
    if (input.peek() == '\n' || input.peek() == ' ') input.get();
    return valid;
}

/**
 * Util function that reads a string from the input stream until a \n is reached
 */
static inline bool getline(std::ifstream &input, std::string &buff) {
    if (std::getline(input, buff, '\n')) return true;
    return false;
}

/**
 * Util function that transforms the input string into his lowered version
 */
static inline void toLower(std::string &s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}


/**
 * Stores .vtk generic data
 */
class DataSection {
public:

    /**
     * Construct a named data field with uninitialized data
     */
    explicit DataSection(std::string section_name) : _section_name(std::move(section_name)) {
        // ensures name is underscore spaced
        std::replace(_section_name.begin(),
                     _section_name.end(),
                     ' ', '_');
    }

    /**
     * Prints .vtk binary formatted content of the class
     */
    virtual bool writeSection(std::ostream &output) = 0;

    /**
     * Reads .vtk binary formatted content of the class
     */
    virtual bool readSection(std::ifstream &input, size_t data_size) = 0;

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
 * Stores .vtk scalars array of values
 */
template<typename Value_Type>
class ScalarsSection : public DataSection {

public:
    /**
     * Data entry value type
     */
    typedef typename std::vector<Value_Type> entry_T;

    /**
      * Construct a named data field with uninitialized data
      */
    explicit ScalarsSection(std::string section_name) : DataSection(section_name) {}

    /**
     * Construct a named data field with initialized data
     */
    ScalarsSection(std::string section_name, std::vector<entry_T> &values)
            : DataSection(section_name) {
        setValues(values);
    }

    // Override functions -------------------------------------------------------- //
    bool writeSection(std::ostream &output) override {
        // standard data field header
        output << "SCALARS " << _section_name << " " << _type_name << " " << _entry_size << "\n";
        output << "LOOKUP_TABLE default\n";

        // binary of all values
        for (const auto &entry: this->_values) {
            for (auto val: entry) {
                swap_endian(val);
                // print binary
                output.write(reinterpret_cast<char *>(&val), sizeof(Value_Type));
            }
        }
        output << "\n";

        return true;
    }

    bool readSection(std::ifstream &input, size_t data_size) override {
        std::string temp;

        getword(input, temp);
        toLower(temp);
        if (temp != "scalars") {
            fprintf(stderr, "DataSection: Trying to read wrong data section type");
            return false;
        }

        getword(input, _section_name);
        getword(input, temp);
        if ((_type_name == "double" || _type_name == "float") && (temp != "double" && temp != "float") ||
            (_type_name == "int" && temp != "int")) {
            fprintf(stderr, "DataSection: Trying yo read wrong data type");
            return false;
        }

        // pick the right reading and converting function based on class type and file-data-type
        std::function<Value_Type()> valueReader;
        if (_type_name != temp) {
            if (temp == "double")
                // in this case file has type double but class is float, it's ok, but it needs to be converted
                valueReader = [&input]() {
                    double val;
                    input.read(reinterpret_cast<char *>(&val), sizeof(double));
                    swap_endian(val);
                    return static_cast<Value_Type>(val);
                };
            else
                // this case is the opposite of the previous
                valueReader = [&input]() {
                    float val;
                    input.read(reinterpret_cast<char *>(&val), sizeof(float));
                    swap_endian(val);
                    return static_cast<Value_Type>(val);
                };
        } else {
            // in this case the type of the file and class are the same, so no conversion
            valueReader = [&input]() {
                Value_Type val;
                input.read(reinterpret_cast<char *>(&val), sizeof(Value_Type));
                swap_endian(val);
                return val;
            };
        }

        getword(input, _entry_size);

        // read default lookuptable
        std::getline(input, temp, '\n');

        // reads all values
        for (size_t i = 0; i < data_size; ++i) {
            auto *t = new std::vector<Value_Type>();
            for (size_t j = 0; j < entry_size; ++j)
                t->push_back(valueReader());
            _values.push_back(*t);
        }

        // skips until next section
        getline(input, temp);

        return true;
    }

    size_t getSize() override { return _values.size(); }
    // ---------------------------------------------------------------------------- //

    /**
     * Setup the section values
     */
    ScalarsSection &setValues(std::vector<entry_T> &values) {

        bool values_coherent = true;

        size_t s = values[0].size();
        for (auto entry: values)
            if (entry.size() != s) {
                values_coherent = false;
                break;
            }

        if (values_coherent)
            _values = values;

        return *this;
    }

    /**
     * Alias for entry size value
     */
    const size_t &entry_size = _entry_size;

    /**
     * Get a copy of values
     */
    std::vector<entry_T> getValues() { return std::move(_values); }

protected:

    /**
     * String representation of the data-field-entry-values type
     */
    const std::string _type_name = getTypeName<Value_Type>();

    /**
     * Data entry list
     */
    std::vector<entry_T> _values;

    /**
     * Size of tuples that each value is composed by
     */
    size_t _entry_size = 0;
};


/**
 * Stores and export .vtk formatted file
 */
class VTKFile {
public:
    /**
     * Build a class with the given file description
     */
    explicit VTKFile(std::string description = "") : _description(std::move(description)) {
        // ensure description maximum length
        if (_description.length() > 255)
            _description = _description.substr(0, 255) + "\n";
    }

    /**
     * Construct a class with the given file description
     * Initialize the dataset with the given physical dimension and node number
     */
    VTKFile(std::array<Real, 3> physical_dimension, std::array<size_t, 3> nodes, std::string description = "")
            : VTKFile(std::move(description)) {
        setupDataset(physical_dimension, nodes);
    }

    /**
     * Initialize the dataset with the given physical dimension and node number
     */
    void setupDataset(std::array<Real, 3> physical_dimension, std::array<size_t, 3> n) {
        _nodes = n;
        _spacing = {physical_dimension[0] / static_cast<Real>(_nodes[0]),
                    physical_dimension[1] / static_cast<Real>(_nodes[1]),
                    physical_dimension[2] / static_cast<Real>(_nodes[2])};
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
                fprintf(stderr, "VTKFile: Tried to insert data that doesn't match dataset\n");
            }
        } else {
            fprintf(stderr, "VTKFile: Tried to insert data before setup dataset\n");
        }
        return *this;
    }

    /**
     * Append a list of data section to the file
     */
    VTKFile &operator<<(const std::vector<DataSection *> &dataSections) {
        for (DataSection *data: dataSections)
            this->operator<<(*data);

        return *this;
    }

    /**
     * Print the content of the class into a file with the given name
     */
    bool writeFile(std::string file_name) {
        // ensures name is underscore spaced
        std::replace(file_name.begin(),
                     file_name.end(),
                     ' ', '_');

        // ensure file extension
        if (!file_name.ends_with(".vtk"))
            file_name += ".vtk";


        if (!_valid_dataset || !_valid_data) {
            fprintf(stderr, "VTKFile: {%s} Tried to print file before setup\n", file_name.c_str());
            return false;
        }

        std::ofstream file;
        file.open(file_name, std::ios::out | std::ios::binary);
        if (!file.is_open()) {
            fprintf(stderr, "VTKFile: {%s} Unable to open the output file\n", file_name.c_str());
            return false;
        }

        writeHeader(file);

        writeDataset(file);

        file << "POINT_DATA " << _data_size << "\n";
        for (DataSection *data: _data_sections)
            data->writeSection(file);


        return true;
    }

    /**
     * Reads a file with the given name and loads class values
     */
    bool readFile(const std::string &file_name) {
        std::ifstream file;
        file.open(file_name, std::ios::in | std::ios::binary);
        if (!file.is_open()) {
            fprintf(stderr, "VTKFile: {%s} Unable to open the input file\n", file_name.c_str());
            return false;
        }

        if (!readHeader(file)) return false;
        if (!readDataset(file)) return false;
        if (!readData(file)) return false;

        return true;
    }


    /**
     * Alias for dataset nodes number
     */
    const std::array<size_t, 3> &nodes = _nodes;

    /**
     * Alias for datased nodes spacing
     */
    const std::array<Real, 3> &spacing = _spacing;

private:
    /**
     * Description of the file
     */
    std::string _description;
    /**
     * Nodes number of dataset
     */
    std::array<size_t, 3> _nodes{};
    /**
     * Node spacing of the dataset
     */
    std::array<Real, 3> _spacing{};
    /**
     * Size of dataset (total number of nodes)
     */
    size_t _data_size{};
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
    std::ostream &writeHeader(std::ostream &output) {
        output << "# vtk DataFile Version 2.0\n";
        output << _description << "\n";
        output << "BINARY\n";
        return output;
    }

    /**
     * Prints the dataset representation into the given output stream
     */
    std::ostream &writeDataset(std::ostream &output) {
        output << "DATASET STRUCTURED_POINTS\n";
        output << "DIMENSIONS " << _nodes[0] << " " << _nodes[1] << " " << _nodes[2] << "\n";
        output << "SPACING " << _spacing[0] << " " << _spacing[1] << " " << _spacing[2] << "\n";
        output << "ORIGIN 0.0 0.0 0.0\n";
        return output;
    }

    /**
     * Reads a .vtk binary formatted header from input stream
     */
    inline bool readHeader(std::ifstream &input) {
        std::string temp;

        // read version
        std::getline(input, temp, '\n');
        toLower(temp);
        // check if version is well-defined
        if (temp != "# vtk datafile version 2.0") {
            fprintf(stderr, "VTKFile: Unable to read header from input file\n");
            return false;
        }

        // read description
        std::getline(input, _description, '\n');

        // read file type
        std::getline(input, temp, '\n');
        toLower(temp);

        // check if type is Binary
        if (temp != "binary") {
            fprintf(stderr, "VTKFile: Unable to read header from input file\n");
            _description = "";
            return false;
        }

        return true;
    }

    /**
     * Reads a .vtk binary formatted dataset from input stream
     */
    inline bool readDataset(std::ifstream &input) {
        std::string temp;

        // read dataset type
        std::getline(input, temp, '\n');
        toLower(temp);
        // check if dataset is structured points
        if (temp != "dataset structured_points") {
            fprintf(stderr, "VTKFile: Unable to read dataset from input file\n");
            return false;
        }

        bool dim = false, sp = false, ori = false;

        // reads all dataset components
        for (int i = 0; i < 3; ++i) {
            input >> temp;
            toLower(temp);
            if (temp == "origin") {
                //TODO add origin values (needed for multiprocessing)
                std::getline(input, temp, '\n');
                ori = true;
            } else if (temp == "dimensions") {
                getword(input, _nodes[0]);
                getword(input, _nodes[1]);
                getword(input, _nodes[2]);
                dim = true;
            } else if (temp == "spacing") {
                getword(input, _spacing[0]);
                getword(input, _spacing[1]);
                getword(input, _spacing[2]);
                sp = true;
            } else {
                fprintf(stderr, "VTKFile: Unable to read dataset from input file\n");
                return false;
            }
        }

        // check if dataset has been correctly read
        if (!sp || !dim || !ori) {
            fprintf(stderr, "VTKFile: Unable to read dataset from input file\n");
            return false;
        }

        _valid_dataset = true;

        return true;
    }

    /**
     * Reads .vtk binary formatted data sections from input stream
     */
    bool readData(std::ifstream &input) {
        std::string temp;

        // check if data is related to points
        getword(input, temp);
        toLower(temp);
        if (temp != "point_data") {
            fprintf(stderr, "VTKFile: Unable to read data from input file\n");
            return false;
        }

        // reads data size
        getword(input, _data_size);

        // reads all data sections until end of file
        while (!input.eof()) {

            //take curr pos
            long pos = input.tellg();

            //read section type
            if (!getword(input,temp)) break;
            toLower(temp);

            //back to pos
            input.seekg(pos, std::ifstream::beg);

            //leave reading to data-section classes
            if (temp == "scalars") {
                DataSection *out = new ScalarsSection<Real>("");

                //tries with scalars of type double/float
                if (out->readSection(input, _data_size)) {
                    this->operator<<(*out);
                } else {
                    //otherwise back to pos and re-tries with scalars of type int
                    input.seekg(pos, std::ifstream::beg);

                    out = new ScalarsSection<int>("");
                    if (out->readSection(input, _data_size)) {
                        this->operator<<(*out);
                    } else {
                        //if no-one works there has been an error
                        return false;
                    }
                }
            } else {
                fprintf(stderr, "VTKFile: Unable to read data from input file\n");
                return false;
            }
        }

        return true;
    }

};


#endif //AEROHPC_A_VTK_H
