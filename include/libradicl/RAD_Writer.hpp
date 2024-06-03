
#ifndef LIBRADICL_RAD_WRITER_HPP
#define LIBRADICL_RAD_WRITER_HPP



#include "Buffer.hpp"
#include "Tags.hpp"
#include "Header.hpp"

#include <cstdint>
#include <cstddef>
#include <string>
#include <fstream>


namespace RAD
{

class Single_End_Read;


class RAD_Writer
{
private:

    const Header header;
    const Tag_Defn tag_defn;

    static constexpr std::size_t buf_cap_default = 4 * 1024 * 1024;    // 4 KB.
    Buffer buf;
    uint32_t read_c_in_buf;

    std::ofstream output;


    void flush_chunk();


public:

    explicit RAD_Writer(const Header& header, const Tag_Defn& tag_defn, const std::string& op_file_path, std::size_t buf_cap = buf_cap_default);

    void add(const Single_End_Read& read_rec);

    void close();
};

}



#endif
