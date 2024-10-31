#ifndef AEROHPC_A_LOGGER_H
#define AEROHPC_A_LOGGER_H

#include <iostream>
#include <string>

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


public:
    explicit Logger(size_t bufsize) : _bufsize(bufsize), sectOpen(false), tableOpen(false){}

    Logger &openSection(std::string title) {
        if (sectOpen) closeSection();

        sectOpen = true;

        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev;

        std::string prevS;
        std::string postS;

        for (int i = 0; i < prev; i++) prevS += "═";
        for (int i = 0; i < next; i++) postS += "═";


        out << "╔" << prevS << " " << title << " " << postS << "╗\n";

        return *this;
    }

    Logger &closeSection() {
        if (!sectOpen) return *this;

        sectOpen = false;

        size_t sectSize = _bufsize - 2;

        std::string prevs;

        for (int i = 0; i < sectSize; i++) prevs += "═";

        out << "╚" << prevs << "╝\n";

        return *this;
    }

    Logger &printTitle(std::string title) {
        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev;

        std::string prevs;
        std::string posts;

        for (int i = 0; i < prev; i++) prevs += "─";
        for (int i = 0; i < next; i++) posts += "─";


        out << "║" << prevs << " " << title << " " << posts << "║\n";
        return *this;
    }

    Logger &printTitle(std::string title, Real value) {
        sanitize(title);

        size_t prev = ((_bufsize - 2) - (title.length() + 2)) / 2;
        size_t next = (_bufsize - 2) - (title.length() + 2) - prev - 10 - 2 - 2;

        std::string prevs;
        std::string posts;

        for (int i = 0; i < prev; i++) prevs += "─";
        for (int i = 0; i < next; i++) posts += "─";


        out << "║" << prevs << " " << title << " " << posts << " ";
        printf("%0.4e", value);
        out << " ──║\n";
        return *this;
    }

    Logger &printValue(int formatRatio, const std::string &value_name, Real value) {

        size_t prev = ((_bufsize - 4) / formatRatio) - value_name.length();
        size_t next = ((_bufsize - 4)) - prev - value_name.length() - 2 - 10;

        std::string prevS;
        std::string nextS;

        for (int i = 0; i < prev; i++) prevS += " ";
        for (int i = 0; i < next; i++) nextS += " ";

        out << "║ " << prevS << value_name << ": ";
        printf("%0.4e", value);
        out << nextS << " ║\n";

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

        std::string prevS;
        std::string nextS;


        for (int i = 0; i < prev; i++) prevS += " ";
        for (int i = 0; i < next; i++) nextS += " ";

        out << "║ " << prevS << value_name << ": " << value << nextS << " ║\n";

        return *this;
    }

    Logger &printValue(int formatRatio, const std::string &value_name, const std::string &str) {

        size_t prev = ((_bufsize - 4) / formatRatio) - value_name.length();
        size_t next = ((_bufsize - 4)) - prev - value_name.length() - 2 - str.length();

        std::string prevS;
        std::string nextS;

        for (int i = 0; i < prev; i++) prevS += " ";
        for (int i = 0; i < next; i++) nextS += " ";

        out << "║ " << prevS << value_name << ": " << str << nextS << " ║\n";

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

        for (int i = 0; i < colSize; i++) endTable += "─";
        for (int i = 0; i < p; i++) tLine += "─";
        tLine += " ";
        tLine += colNames[0];
        tLine += " ";
        for (int i = 0; i < n; i++) tLine += "─";
        nchar += p + 2 + colNames[0].size() + n;

        for (int j = 1; j < colNames.size(); ++j) {
            endTable += "┴";
            tLine += "┬";

            p = (colSize - 2 - colNames[j].size()) / 2;
            n = colSize - 2 - colNames[j].size() - p;

            for (int i = 0; i < colSize; i++) endTable += "─";
            for (int i = 0; i < p; i++) tLine += "─";
            tLine += " ";
            tLine += colNames[j];
            tLine += " ";
            for (int i = 0; i < n; i++) tLine += "─";
            nchar += 1 + p + 2 + colNames[j].size() + n;
        }

        tLine += "┐";
        endTable += "┘";

        nchar++;

        size_t prev = ((_bufsize - 2) - nchar) / 2;
        size_t next = (_bufsize - 2) - nchar - prev;

        std::string prevs;
        std::string posts;

        for (int i = 0; i < prev; i++) prevs += " ";
        for (int i = 0; i < next; i++) posts += " ";


        endTable = "║" + prevs + endTable + posts + "║\n";
        out << "║" << prevs << tLine << posts << "║\n";

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

        size_t colSize = ((_bufsize -4 - 1) / (values.size() + 1)) - 1;

        std::string tLine = "│";
        size_t nchar = 1;

        size_t p = (colSize - 2 - idSize) / 2;
        size_t n = colSize - 2 - idSize - p;
        for (int i = 0; i < p; i++) tLine += " ";
        tLine += " ";
        tLine += to_string(rowId);
        tLine += " ";
        for (int i = 0; i < n; i++) tLine += " ";
        nchar += p + 2 + idSize + n;


        nchar += 1 + values.size() * (1 + colSize);

        size_t prev = ((_bufsize - 2) - nchar) / 2;
        size_t next = (_bufsize - 2) - nchar - prev;

        std::string prevs = "║";
        std::string posts;

        for (int i = 0; i < prev; i++) prevs += " ";
        for (int i = 0; i < next; i++) posts += " ";

        out << prevs << tLine;
        tLine = "";

        for (double value: values) {
            tLine += "│";
            p = (colSize - 2 - 10) / 2;
            n = colSize - 2 - 10 - p;
            for (int i = 0; i < p; i++) tLine += " ";
            tLine += " ";
            out << tLine;
            printf("%0.4e", value);
            tLine = " ";
            for (int i = 0; i < n; i++) tLine += " ";
        }

        tLine += "│";

        out << tLine << posts << "║\n";

        return *this;
    }

    Logger &empty() {
        out << endl;
        return *this;
    }
};

#endif //AEROHPC_A_LOGGER_H
