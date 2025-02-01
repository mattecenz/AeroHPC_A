#ifndef AEROHPC_A_LOGGER_HPP
#define AEROHPC_A_LOGGER_HPP

#include <iostream>
#include <string>
#include <iterator>
#include <sstream>
#include <iomanip>


/**
* Macro for retrieving symmetric spacing
*/
#define getPrevNext(name, strlen, size) \
const size_t name##prev = (size - strlen) / 2;  \
const size_t name##next = size - strlen - name##prev

class Logger {
    /**
    * Symmetric padding size
    */
    size_t paddingS = 2;

    /**
    * Output stream of logger
    */
    std::ostream &out = std::cout;

    /**
    * Horizontal size of the logger
    */
    size_t _buffSize;

    /**
    * If the logger has a section open
    */
    bool sectOpen;

    /**
    * If the logger has a table open
    */
    bool tableOpen;

    /**
    * Buffered footer of a table
    */
    std::string endTable;

    /**
    * Check whether the input string is shorter than the given length,
    * otherwise transforms it into a dotted format
    */
    static void sanitize(std::string string, const size_t length) {
        if (string.length() > length) {
            string.resize(length - 1);
            string += ".";
        }
    }

    /**
    * Returns a string composed by repeating num times the input string
    * If the returned string is longer than the logger size then returns an empty string
    */
    std::string repeat(const std::string &input, const size_t num) const {
        const size_t _num = num > _buffSize ? 0 : num;
        std::stringstream os;
        std::fill_n(std::ostream_iterator<std::string>(os), _num, input);
        return os.str();
    }

public:
    /**
    * Construct a logger of at least 80 characters, with no section nor table opened
    */
    explicit Logger(size_t bufsize) : sectOpen(false), tableOpen(false) {
        if (bufsize < 80) _buffSize = 80;
        else _buffSize = bufsize;
    }

    /**
    * Open a section and print a section header line with given title
    */
    Logger &openSection(const std::string &title) {
        if (sectOpen) closeSection();

        sectOpen = true;

        sanitize(title, _buffSize - paddingS - paddingS);

        getPrevNext(sect, title.size(), _buffSize - paddingS - paddingS);

        out << "╔" << repeat("═", sectprev) << " " << title << " " << repeat("═", sectnext) << "╗\n";

        return *this;
    }

    /**
    * Close a section and print a section footer line
    */
    Logger &closeSection() {
        if (!sectOpen) return *this;

        sectOpen = false;

        size_t sectSize = _buffSize - 2;

        out << "╚" << repeat("═", sectSize) << "╝\n";

        return *this;
    }

    /**
    * Print a title line with the given title
    */
    Logger &printTitle(const std::string &title) {
        sanitize(title, _buffSize - paddingS - paddingS);

        getPrevNext(t, title.size(), _buffSize - paddingS - paddingS);

        out << "║" << repeat("─", tprev) << " " << title << " " << repeat("─", tnext) << "║\n";

        return *this;
    }

    /**
    * Print a title line with the given title and value
    */
    Logger &printTitle(const std::string &title, const Real value, const long int precision = 4) {
        //                        0.  precision  e+xx  spaces     lines
        const size_t valueSpace = 2 + precision + 4 + paddingS + paddingS;
        const size_t titleSpace = (((_buffSize - paddingS - paddingS) / 2) - valueSpace) * 2;

        if (titleSpace <= 0) return *this;

        sanitize(title, titleSpace);

        getPrevNext(t, title.size(), _buffSize - paddingS - paddingS);

        out << "║" << repeat("─", tprev) << " " << title << " " << repeat("─", tnext - valueSpace + 1) << " "
                << std::setprecision(precision) << std::scientific << value << " ─║\n";

        return *this;
    }

    /**
    * Print a value information,
    * formatRatio is the portion of line the info will be placed on,
    * valueName is the information name,
    * value is the information value,
    * precision is the decimal size of the value
    */
    Logger &printValue(const long int formatRatio, const std::string &valueName, const Real value, const long int precision = 4) {
        const size_t valueSpace = 2 + precision + 4 + paddingS;
        const size_t titleSpace = (_buffSize - paddingS) / formatRatio - paddingS - 1;

        if (titleSpace <= 0) return *this;
        if ((_buffSize - paddingS) * (1 - 1 / formatRatio) < valueSpace) return *this;

        sanitize(valueName, titleSpace);

        getPrevNext(t, valueName.size(), (_buffSize - paddingS) / formatRatio);
        const size_t vNext = (_buffSize - paddingS) - ((_buffSize - paddingS) / formatRatio) - valueSpace;

        out << "║" << repeat(" ", tprev + tnext - 2) << valueName << ":  "
                << std::setprecision(precision) << std::scientific << value << repeat(" ", vNext) << " ║\n";

        return *this;
    }

    /**
    * Print a value information,
    * formatRatio is the portion of line the info will be placed on,
    * valueName is the information name,
    * value is the information value
    */
    Logger &printValue(const long int formatRatio, const std::string &valueName, const long int value) {
        long int copy = value;
        size_t valueSize = 1;
        while (copy / 10 > 0) {
            valueSize++;
            copy /= 10;
        }

        const size_t valueSpace = valueSize + paddingS;
        const size_t titleSpace = (_buffSize - paddingS) / formatRatio - paddingS - 1;

        if (titleSpace <= 0) return *this;
        if ((_buffSize - paddingS) * (1 - 1 / formatRatio) < valueSpace) return *this;

        sanitize(valueName, titleSpace);

        getPrevNext(t, valueName.size(), (_buffSize - paddingS) / formatRatio);
        const size_t vNext = (_buffSize - paddingS) - ((_buffSize - paddingS) / formatRatio) - valueSpace;

        out << "║" << repeat(" ", tprev + tnext - 2) << valueName << ":  "
                << std::dec << value << repeat(" ", vNext) << " ║\n";

        return *this;
    }

    /**
    * Print a value information,
    * formatRatio is the portion of line the info will be placed on,
    * valueName is the information name,
    * value is the information value
    */
    Logger &printValue(const long int formatRatio, const std::string &valueName, const int value) {
        long int copy = value;
        size_t valueSize = 1;
        while (copy / 10 > 0) {
            valueSize++;
            copy /= 10;
        }

        const size_t valueSpace = valueSize + paddingS;
        const size_t titleSpace = (_buffSize - paddingS) / formatRatio - paddingS - 1;

        if (titleSpace <= 0) return *this;
        if ((_buffSize - paddingS) * (1 - 1 / formatRatio) < valueSpace) return *this;

        sanitize(valueName, titleSpace);

        getPrevNext(t, valueName.size(), (_buffSize - paddingS) / formatRatio);
        const size_t vNext = (_buffSize - paddingS) - ((_buffSize - paddingS) / formatRatio) - valueSpace;

        out << "║" << repeat(" ", tprev + tnext - 2) << valueName << ":  "
                << std::dec << value << repeat(" ", vNext) << " ║\n";

        return *this;
    }

    /**
    * Print a value information,
    * formatRatio is the portion of line the info will be placed on,
    * valueName is the information name,
    * str is the information string description
    */
    Logger &printValue(const int formatRatio, const std::string &valueName, const std::string &str) {
        const size_t valueSpace = str.size() + paddingS;
        const size_t titleSpace = (_buffSize - paddingS) / formatRatio - paddingS - 1;

        if (titleSpace <= 0) return *this;
        if ((_buffSize - paddingS) * (1 - 1 / formatRatio) < valueSpace) return *this;

        sanitize(valueName, titleSpace);

        getPrevNext(t, valueName.size(), (_buffSize - paddingS) / formatRatio);
        const size_t vNext = (_buffSize - paddingS) - ((_buffSize - paddingS) / formatRatio) - valueSpace;

        out << "║" << repeat(" ", tprev + tnext - 2) << valueName << ":  "
                << str << repeat(" ", vNext) << " ║\n";

        return *this;
    }

    /**
    * Print a value information,
    * formatRatio is the portion of line the info will be placed on,
    * valueName is the information name,
    * value is the information value
    */
    Logger &printValue(const long int formatRatio, const std::string &valueName, const bool value) {
        const std::string value_str = value ? "true" : "false";
        const size_t valueSpace = value_str.size() + paddingS;
        const size_t titleSpace = (_buffSize - paddingS) / formatRatio - paddingS - 1;

        if (titleSpace <= 0) return *this;
        if ((_buffSize - paddingS) * (1 - 1 / formatRatio) < valueSpace) return *this;

        sanitize(valueName, titleSpace);

        getPrevNext(t, valueName.size(), (_buffSize - paddingS) / formatRatio);
        const size_t vNext = (_buffSize - paddingS) - ((_buffSize - paddingS) / formatRatio) - valueSpace;

        out << "║" << repeat(" ", tprev + tnext - 2) << valueName << ":  "
                << value_str << repeat(" ", vNext) << " ║\n";

        return *this;
    }

    /**
    * Open a table and print a table header given the column names,
    * rowIdName is the name of the first column,
    * colNames is a collection of column names
    */
    Logger &openTable(const std::string &rowIdName, const std::vector<std::string> &colNames) {
        if (tableOpen) closeTable();

        const size_t colSize = (_buffSize - paddingS - paddingS - 1) / (colNames.size() + 1)  - 1;

        std::stringstream tOpen;
        std::stringstream tClose;
        tOpen << "┌";
        tClose << "└";

        sanitize(rowIdName, colSize - paddingS);
        getPrevNext(n0, rowIdName.size(), colSize - paddingS);

        tOpen << repeat("─", n0prev) << " " << rowIdName << " " << repeat("─", n0next);
        tClose << repeat("─", n0prev + 1 + rowIdName.size() + 1 + n0next);

        for (int j = 0; j < colNames.size(); ++j) {
            sanitize(colNames[j], colSize - paddingS);
            getPrevNext(nj, colNames[j].size(), colSize - paddingS);

            tOpen << "┬" << repeat("─", njprev) << " " << colNames[j] << " " << repeat("─", njnext);
            tClose << "┴" << repeat("─", njprev + 1 + colNames[j].size() + 1 + njnext);
        }
        tOpen << "┐";
        tClose << "┘";

        const size_t tableSize = (colSize + 1) * (colNames.size() + 1) + 1;
        getPrevNext(t, tableSize, _buffSize - paddingS);

        out << "║" << repeat(" ", tprev) << tOpen.str() << repeat(" ", tnext) + "║\n";

        std::stringstream endTableS;
        endTableS << "║" << repeat(" ", tprev) << tClose.str() << repeat(" ", tnext) + "║\n";
        endTable = endTableS.str();

        tableOpen = true;
        return *this;
    }

    /**
    * Close a table and print the table footer
    */
    Logger &closeTable() {
        if (!tableOpen) return *this;

        out << endTable;

        tableOpen = false;

        return *this;
    }

    /**
    * Print a table line,
    * rowId is the first column value,
    * values is a collection of values,
    * precision is the decimal size of values
    */
    Logger &printTableValues(const long int rowId, const std::vector<Real> &values, const long int precision = 4) {
        long int copyId = rowId;
        size_t idSize = 1;
        while (copyId / 10 > 0) {
            idSize++;
            copyId /= 10;
        }

        const size_t colSize = (_buffSize - paddingS - paddingS - 1) / (values.size() + 1) - 1;

        size_t idSpace = idSize + paddingS;
        const size_t valueSpace = 2 + precision + 4 + paddingS;

        if (idSpace > colSize) {
            copyId = -1;
            idSpace = 4;
        } else copyId = rowId;

        size_t tableSize = (colSize + 1) * (values.size() + 1) + 1;

        std::stringstream tLine;
        tLine << "│";

        getPrevNext(id, idSpace, colSize);
        getPrevNext(v, valueSpace, colSize);


        tLine << repeat(" ", idprev + 1) << copyId << repeat(" ", idnext + 1);

        for (int j = 0; j < values.size(); ++j) {
            tLine << "│" << repeat(" ", vprev + 1) << std::setprecision(precision) << std::scientific << values[j] << repeat(" ", vnext + 1);
        }
        tLine << "│";

        getPrevNext(t, tableSize, _buffSize - paddingS);

        out << "║" << repeat(" ", tprev) << tLine.str() << repeat(" ", tnext) + "║\n";

        return *this;
    }

    /**
    * Print an empty line into a section
    */
    Logger &spacer() {
        out << "║" << repeat(" ", _buffSize - 2) << "║\n";
        return *this;
    }

    /**
    * Print an empty line
    */
    Logger &empty() {
        out << std::endl;
        return *this;
    }
};

/**
 * Define global logger
 */
inline Logger logger(150);

#define enabledLogger if(IS_MAIN_PROC) logger

#endif //AEROHPC_A_LOGGER_HPP
