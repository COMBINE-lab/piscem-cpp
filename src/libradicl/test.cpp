
#include "RAD_Writer.hpp"
#include "Alignment_Record.hpp"
#include "Read_Record.hpp"
#include "Byte_Array.hpp"
#include "Tags.hpp"
#include "reference_index.hpp"


int main()
{
    const std::vector<std::string> refs{"ref-path-0"};
    const RAD::Header header(false, refs.size(), refs);

    RAD::Tag_Defn tag_defn;
    tag_defn.add_file_tag<RAD::Type::str>("read_tech");

    tag_defn.add_read_tag<RAD::Type::str>("read_name");
    tag_defn.add_read_tag<RAD::Type::f32>("avg_read_qual");

    tag_defn.add_aln_tag<RAD::Type::u32>("aln_pos");
    tag_defn.add_aln_tag<RAD::Type::f32>("aln_score");

    RAD::RAD_Writer rad_writer(header, tag_defn, "op-RAD-path");

    RAD::Single_End_Read read_rec;
    RAD::Aln_Record aln_rec;
    constexpr std::size_t read_example = 10;
    constexpr std::size_t aln_count = 5;
    for(std::size_t i = 0; i < read_example; ++i)
    {
        read_rec.set(aln_count, 100);
        read_rec.add_tag(RAD::Type::str("dummy-read-name"));
        read_rec.add_tag(RAD::Type::f32(20.5));

        for(std::size_t j = 0; j < aln_count; ++j)
        {
            aln_rec.set(0, 1);
            aln_rec.add_tag(RAD::Type::u32(77));
            aln_rec.add_tag(RAD::Type::f32(99.0));

            read_rec.add_aln_rec(aln_rec);
        }

        rad_writer.add(read_rec);
    }


    rad_writer.close();

    return 0;
}
