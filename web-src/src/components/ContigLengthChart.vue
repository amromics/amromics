<style scoped>
</style>
<template>
<div id="contig_chart">
</div>
</template>
<script>
/* eslint-disable */
import {N50Chart} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'ContigLengthChart',
    props: ["list_contig"],
    data() {
        return {
            loading: false,

        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('contig_chart')
      //console.log(this.list_contig);
      var n50chart = new N50Chart(ctx);
      n50chart.load(this.list_contig);
      n50chart.setOptions({width:1200,height:550});
      n50chart.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
