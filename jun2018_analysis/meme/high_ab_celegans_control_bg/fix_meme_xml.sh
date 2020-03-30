#!/usr/bin/env bash
MEMEPREFIX="/home/hisakatha/memesuite_5_1_1_install"
cd result/meme_out &&
    cp -f meme.xml meme.back.xml &&
    if [[ -e meme.html ]]; then cp -f meme.html meme.back.html; fi &&
    sed -E -e "s@([ \t]*<maxsites>).+(</maxsites>[ \t]*)@\10\2@" meme.back.xml > meme.xml &&
    $MEMEPREFIX/libexec/meme-5.1.1/meme_xml_to_html meme.xml meme.html
