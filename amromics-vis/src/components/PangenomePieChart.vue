<style scoped>
</style>
<template>
<div id="piechart">
</div>
</template>
<script>
/* eslint-disable */
import {PanPieChart} from "@/amromicsjs";
// import SampleIGV from "@/components/Visualization/IGV";
export default {
    name: 'PangenomePieChart',
    props: ['core_data'],
    data() {
        return {
            loading: false,
          
        };
    },
    mounted() {
      this.loading = true;
      var ctx=document.getElementById('piechart');
    console.log(this.core_data);
      var pie = new PanPieChart(ctx);
      pie.load(this.core_data);
      
      pie.draw();
      this.loading = false;
    },
    async created() {

      
    },
    method: {
      
    }
};
</script>
