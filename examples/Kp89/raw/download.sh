#!/bin/bash

for ac in ERR349747 ERR349756 ERR349757 ERR349758 ERR349759 ERR349763 ERR349764 ERR349765 ERR349748 ERR349766 ERR349767 ERR349768 ERR349769 ERR349770 ERR349773 ERR349774 ERR349775 ERR349749 ERR349777 ERR349778 ERR349780 ERR349784 ERR349785 ERR349786 ERR349787 ERR349788 ERR349789 ERR349790 ERR349791 ERR349795 ERR349751 ERR349797 ERR349798 ERR349799 ERR349800 ERR349803 ERR349804 ERR349805 ERR349752 ERR349806 ERR349807 ERR349808 ERR349811 ERR349812 ERR349813 ERR349814 ERR349815 ERR349753 ERR349816 ERR349817 ERR349818 ERR349819 ERR349820 ERR349821 ERR349822 ERR349823 ERR349824 ERR349825 ERR349754 ERR349826 ERR349828 ERR349831 ERR349833 ERR349834 ERR349835 ERR349755 ERR349838 ERR349839 ERR349841 ERR349842 ERR349843 ERR349844 ERR349845 ERR349846 ERR349847 ERR349848 ERR349849 ERR349852 ERR349853 ERR349854 ERR317538 ERR317539 ERR317529 ERR317530 ERR317531 ERR317532 ERR317533 ERR317535 ERR317536;do
    if [ ! -f ${ac}_2.fastq.gz ];then
        echo " Downloading Accession ${ac} ..."
        fasterq-dump --progress --split-3 ${ac} && gzip ${ac}_1.fastq && gzip ${ac}_2.fastq
    else
        echo " Accession ${ac} has been downloaded!"
    fi
done




