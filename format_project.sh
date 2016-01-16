find ./ \( -iname "*.hpp" -o -iname "*.cpp" -o -iname "*.cu" \) | xargs clang-format -i
