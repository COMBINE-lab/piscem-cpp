
#include "../../include/libradicl/RAD_Writer.hpp"
#include "../../include/libradicl/Read_Record.hpp"


namespace RAD
{

RAD_Writer::RAD_Writer(const Header& header, const Tag_Defn& tag_defn, const std::string& op_file_path, const std::size_t buf_cap):
      header(header)
    , tag_defn(tag_defn)
    , buf(buf_cap, output)
    , read_c_in_buf(0)
    , output(op_file_path, std::ios::out | std::ios::binary)
{
    header.write(buf);
    tag_defn.write(buf);
}


void RAD_Writer::add(const Single_End_Read& read_rec)
{
    if(buf.size() + read_rec.size_in_bytes() > buf.capacity())
        flush_chunk();

    assert(buf.size() + read_rec.size_in_bytes() <= buf.capacity());
    buf.add(read_rec.byte_array());
    read_c_in_buf++;
}


void RAD_Writer::flush_chunk()
{
    const Type::u32 chunk_sz = buf.size();
    output.write(reinterpret_cast<const char*>(&chunk_sz), sizeof(chunk_sz));
    output.write(reinterpret_cast<const char*>(&read_c_in_buf), sizeof(read_c_in_buf));

    buf.flush();
    read_c_in_buf = 0;
}


void RAD_Writer::close()
{
    if(!buf.empty())
    {
        assert(read_c_in_buf > 0);
        flush_chunk();
    }
}

}
