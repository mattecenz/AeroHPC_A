#ifndef AEROHPC_A_LOGGER_HPP
#define AEROHPC_A_LOGGER_HPP

#include <iostream>
#include <string>
#include <iterator>
#include <sstring>


class Logger {
private:
    std::ostream &out = std::cout;
    size_t _bufsize;
    bool sectOpen;
    bool tableOpen;
    std::string endTable;

    void sanitize(std::string &string) const {
        if (string.length() + 2 > (_bufsize - 2))
            string.resize(_bufsize - 2);
    }

    static std::string repeat(const std::string &input, size_t num) {
        std::stringstream os;
        std::fill_n(std::ostream_iterator<std::string>(os), num, input);
        return os.str();
    }

    static std::string format(Real &value){
        std::stringstream os;
        os.precision(4);
        os << std::scientific << value;
        return os.str();
    }


public:
    explicit Logger(size_t bufsize) : _bufsize(bufsize), sectOpen(false), tableOpen(false) {}

    Logger &openSection(std::string title) {
        if (sectOpen) closeSection();

        sectOpen = true;

        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev;

        out << "╔" << repeat("═", prev) << " " << title << " " << repeat("═", next) << "╗\n";

        return *this;
    }

    Logger &closeSection() {
        if (!sectOpen) return *this;

        sectOpen = false;

        size_t sectSize = _bufsize - 2;

        out << "╚" << repeat("═", sectSize) << "╝\n";

        return *this;
    }

    Logger &printTitle(std::string title) {
        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev;

        out << "║" << repeat("─", prev) << " " << title << " " << repeat("─", next) << "║\n";
        return *this;
    }

    Logger &printTitle(std::string title, Real value) {
        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev - 10 - 2 - 2;

        out << "║" << repeat("─", prev) << " " << title << " " << repeat("─", next) << " ";
        printf("%0.4e", value);
        out << " ──║\n";
        return *this;
    }

    Logger &printValue(int formatRatio, const std::string &value_name, Real value) {

        size_t prev = ((_bufsize - 4) / formatRatio) - value_name.length();
        size_t next = ((_bufsize - 4)) - prev - value_name.length() - 2 - 10;

        out << "║ " << repeat(" ", prev) << value_name << ": ";
        printf("%0.4e", value);
        out << repeat(" ", next) << " ║\n";

        return *this;
    }

    Logger &printValue(int formatRatio, const std::string &value_name, long int value) {

        long int copy = value;
        size_t valueSize = 1;
        while (copy / 10 > 0) {
            valueSize++;
            copy /= 10;
        }

        size_t prev = ((_bufsize - 4) / formatRatio) - value_name.length();
        size_t next = ((_bufsize - 4)) - prev - value_name.length() - 2 - valueSize;

        out << "║ " << repeat(" ", prev) << value_name << ": " << value << repeat(" ", next) << " ║\n";

        return *this;
    }

    Logger &printValue(int formatRatio, const std::string &value_name, const std::string &str) {

        size_t prev = ((_bufsize - 4) / formatRatio) - value_name.length();
        size_t next = ((_bufsize - 4)) - prev - value_name.length() - 2 - str.length();

        out << "║ " << repeat(" ", prev) << value_name << ": " << str << repeat(" ", next) << " ║\n";

        return *this;
    }

    Logger &openTable(const std::vector<std::string> &colNames) {
        if (tableOpen) closeTable();

        tableOpen = true;

        size_t colSize = ((_bufsize - 4 - 1) / colNames.size()) - 1;

        std::string tLine = "┌";
        endTable = "└";
        size_t nchar = 1;

        size_t p = (colSize - 2 - colNames[0].size()) / 2;
        size_t n = colSize - 2 - colNames[0].size() - p;

        endTable += repeat("─", colSize);
        tLine += repeat("─", p) + " " + colNames[0] + " " + repeat("─", n);
        nchar += p + 2 + colNames[0].size() + n;

        for (int j = 1; j < colNames.size(); ++j) {
            p = (colSize - 2 - colNames[j].size()) / 2;
            n = colSize - 2 - colNames[j].size() - p;

            endTable += "┴" + repeat("─", colSize);
            tLine += "┬" + repeat("─", p) + " " + colNames[j] + " " + repeat("─", n);
            nchar += 1 + p + 2 + colNames[j].size() + n;
        }

        tLine += "┐";
        endTable += "┘";

        nchar++;

        size_t prev = ((_bufsize - 2) - nchar) / 2;
        size_t next = (_bufsize - 2) - nchar - prev;

        endTable = "║" + repeat(" ", prev) + endTable + repeat(" ", next) + "║\n";
        out << "║" << repeat(" ", prev) << tLine << repeat(" ", next) << "║\n";

        return *this;
    }

    Logger &closeTable() {
        if (!tableOpen) return *this;

        tableOpen = false;

        out << endTable;

        return *this;
    }

    Logger &printTableValues(long int rowId, const std::vector<Real> &values) {

        long int copy = rowId;
        size_t idSize = 1;
        while (copy / 10 > 0) {
            idSize++;
            copy /= 10;
        }

        size_t colSize = ((_bufsize - 4 - 1) / (values.size() + 1)) - 1;

        size_t nchar = 1;

        size_t p = (colSize - 2 - idSize) / 2;
        size_t n = colSize - 2 - idSize - p;

        std::string tLine = "│" + repeat(" ", p) + " " + to_string(rowId) + " " + repeat(" ", n);
        nchar += p + 2 + idSize + n;

        for (double value: values) {
            p = (colSize - 2 - 10) / 2;
            n = colSize - 2 - 10 - p;

            tLine += "│" + repeat(" ", p + 1) + format(value) + repeat(" ", n + 1);
        }
        tLine += "│";

        nchar += 1 + values.size() * (1 + colSize);

        size_t prev = ((_bufsize - 2) - nchar) / 2;
        size_t next = (_bufsize - 2) - nchar - prev;

        out << "║" << repeat(" ", prev) << tLine << repeat(" ", next) << "║\n";

        return *this;
    }

    Logger &empty() {
        out << endl;
        return *this;
    }
};

#endif //AEROHPC_A_LOGGER_HPP
